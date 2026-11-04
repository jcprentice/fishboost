{
    library(data.table)
    library(stringr)
    library(purrr)
    library(coda)
    library(HDInterval)
    library(pipebind)

    safe_source("rename_pars.R",
                "utils.R")
}


rebuild_bici_posteriors <- function(dataset = "fb-test",
                                    name = "scen-2-1") {

    if (FALSE) {
        dataset <- "fb-test"
        name <- "scen-2-1"
        dataset <- params$dataset
        name <- params$name
    }

    base_dir <- str_glue("datasets/{dataset}")
    data_dir <- str_glue("{base_dir}/data")
    res_dir  <- str_glue("{base_dir}/results")
    out_dir  <- str_glue("{data_dir}/{name}-out")

    f <- str_glue("{res_dir}/{name}.rds")
    if (file.exists(f)) {
        params <- readRDS(f)$params
        nsample  <- params$nsample
        thinto   <- params$thinto |> min(nsample)
        burnprop <- params$burnprop
    } else {
        # All this to get nsamples and burnin, which are probably just 1e4 and 0.2
        inf_line <- readLines(str_glue("{data_dir}/{name}.bici")) |>
            str_subset("^inference") |>
            str_replace_all(" = ", "=") |>
            str_split_1(" ") |>
            str_subset(".+") |>
            str_split("=")

        inf_pars <- map(inf_line, 2) |> setNames(map(inf_line, 1))

        nsample  <- (inf_pars[["update"]] %||% 1e6L) |> as.integer()
        thinto   <- (inf_pars[["param-output"]] %||% 1e4L) |> as.integer() |> min(nsample)
        burnprop <- (inf_pars[["burn-prop"]] %||% 0.2) |>  as.numeric()
    }

    x <- list.files(str_glue("{out_dir}/output-inf"),
                    "^param_",
                    full.names = TRUE) |>
        str_sort(numeric = TRUE) |>
        map(fread) |>
        map(~ .x[, .SD, .SDcols = !patterns("L\\^|N\\^|Prior")]) |>
        map(~ setnames(.x, rename_bici_pars))

    walk(x, \(xi) {
        if (FALSE) {
            xi <- copy(x[[1]])
        }

        traits <- names(xi) |> str_subset("cov_G") |> str_remove("cov_G_") |> str_1st()

        covs <- str_c(traits, traits)
        cors <- expand.grid(traits, traits) |> apply(1, str_flatten) |>
            matrix(length(traits)) |> bind(x, t(x)[lower.tri(x)])

        walk(covs, \(t1) {
            gt <- str_c("cov_G_", t1)
            et <- str_c("cov_E_", t1)
            pt <- str_c("cov_P_", t1)
            ht <- str_c("h2_", t1)
            xi[, (pt) := get(gt) + get(et)]
            xi[, (ht) := get(gt) / get(pt)]
        })

        walk(cors, \(st) {
            s <- str_sub(st, 1, 1)
            t <- str_sub(st, 2, 2)

            rgst <- str_c("r_G_", st)
            rest <- str_c("r_E_", st)
            rpst <- str_c("r_P_", st)
            cgs <- str_c("cov_G_", s, s)
            cgt <- str_c("cov_G_", t, t)
            cgst <- str_c("cov_G_", st)
            ces <- str_c("cov_E_", s, s)
            cet <- str_c("cov_E_", t, t)
            cest <- str_c("cov_E_", st)
            cps <- str_c("cov_P_", s, s)
            cpt <- str_c("cov_P_", t, t)
            cpst <- str_c("cov_P_", st)

            xi[, (cgst) := get(rgst) * sqrt(get(cgs) * get(cgt))]
            xi[, (cest) := get(rgst) * sqrt(get(ces) * get(cet))]
            xi[, (cpst) := get(cgst) + get(cest)]
            xi[, (rpst) := get(cpst) / sqrt(get(cps) + get(cpt))]
            xi[, c(cgst, cest, cpst) := NULL]
        })


        pnames <- c(str_subset(names(xi), "_P_"),
                    str_subset(names(xi), "h2_"))
        if (length(pnames) > 0) {
            setcolorder(xi, pnames,
                        after = last(str_subset(names(xi), "_[GE]_")),
                        skip_absent = TRUE)
        }
    })

    xs <- x |> map(~ .x[-seq_len(burnprop * .N)]) |> rbindlist()
    xs[, state := .I]

    # Write to trace_combine.tsv
    tc <- str_glue("{out_dir}/trace_combine.tsv")
    message(str_glue("- Writing '{tc}'"))
    fwrite(xs, file = tc, sep = "\t")

    xs[, state := NULL]

    xp <- xs[, .(parameter = .SD |> names(),
                 true_val  = NA_real_,
                 type      = NA_character_,
                 val1      = NA_real_,
                 val2      = NA_real_,
                 mean      = .SD |> map_dbl(mean),
                 median    = .SD |> map_dbl(median),
                 sd        = .SD |> map_dbl(sd),
                 ci95min   = .SD |> map_dbl(quantile, 0.025),
                 ci95max   = .SD |> map_dbl(quantile, 0.975),
                 hdi95min  = .SD |> map(hdi) |> map_dbl("lower"),
                 hdi95max  = .SD |> map(hdi) |> map_dbl("upper"),
                 prior_var = NA_real_,
                 contract  = NA_real_)]

    walk(seq_len(nrow(xp)), \(i) {
        par <- xp$parameter[[i]]
        if (par %in% params$priors$parameter) {
            priors <- params$priors[parameter == par,
                                    .(true_val, type, val1, val2)] |>
                as.list()
            xp[i, c("true_val", "type", "val1", "val2") := priors]
            # set(xp, i, c("true_val", "type", "val1", "val2"),
            #     params$priors[parameter == par, .(true_val, type, val1, val2)] |> as.list())
        }
        if (str_detect(dataset, "sim")) {
            xp[i, true_val := NA_real_]
        }
    })

    inv_var <- function(a, b) {
        (b^2 - a^2) / (2 * log(b /a )) -
            ((b - a) / log(b / a))^2
    }

    xp[type == "uniform", prior_var := (val2 - val1)^2 / 12]
    xp[type == "inverse", prior_var := inv_var(val1, val2)]
    xp[str_starts(parameter, "cov_"), prior_var := 4 * (1 - 2 / pi)]
    xp[str_starts(parameter, "r_"), prior_var := 1 / (2 * 1.2 + 1)]

    xp[, contract := 1 - sd^2 / prior_var]
    xp[, prior_var := NULL]

    if (str_detect(dataset, "sim")) {
        tvs <- with(params$priors, setNames(true_val, parameter))
        xp[parameter %in% names(tvs), true_val := tvs[parameter]]
    }


    # Add in ESS and GR

    # Need to avoid cols with 0 variance
    check <- xs[, map_dbl(.SD, sd)]
    drop_pars <- names(check[check == 0])
    cols <- setdiff(xp$parameter, drop_pars) |>
        str_subset("state|^G_|_P_|h2", negate = TRUE)

    # In case we didn't finish or have an uneven number of samples
    min_rows <- map_int(x, nrow) |> min()

    tfs <- x |>
        map(~ .x[seq_len(min_rows), ..cols]) |>
        map(~ mcmc(.x, start = nsample * burnprop, thin = nsample / thinto))

    grd <- as.mcmc.list(tfs) |> gelman.diag() |> _$psrf[, 2]
    ess <- effectiveSize(tfs) |> round()

    xp[parameter %in% cols, `:=`(ESS = ess, GR = grd)]
    xp[, convergence := fcase(ESS >= 500 & GR < 1.05, "***",
                              ESS >= 200 & GR < 1.1, "**",
                              ESS >= 100 & GR < 1.2, "*",
                              default = "")]

    res_file <- str_glue("{res_dir}/{name}.rds")
    if (file.exists(res_file)) {
        tmp <- readRDS(res_file)
        tmp$parameter_estimates <- xp
        saveRDS(tmp, res_file)
    }

    posteriors_tsv <- str_glue("{out_dir}/posterior.tsv")
    message(str_glue("- Writing '{posteriors_tsv}'"))
    fwrite(xp, posteriors_tsv, sep = "\t")
    xp
}

rebuild_sire_posteriors <- function(dataset = "fb-final",
                                    name = "scen-1-1") {

    if (FALSE) {
        dataset <- "fb-final"
        name <- "scen-1-1"
    }
    # dataset <- params$dataset; name <- params$name

    data_dir <- str_glue("datasets/{dataset}/data")
    out_dir  <- str_glue("{data_dir}/{name}-out")
    res_dir  <- str_glue("datasets/{dataset}/results")

    tc_txt <- str_glue("{out_dir}/trace_combine.txt")
    tc     <- str_glue("{out_dir}/trace_combine.tsv")

    if (file.exists(tc_txt)) file.rename(tc_txt, tc)

    pf <- str_glue("{res_dir}/{name}.rds")
    if (file.exists(pf)) {
        params <- readRDS(pf)$params
        nsample <- params$nsample
        burnprop <- params$burnprop
        thinto <- params$thinto
        priors <- params$priors
    } else {
        lines <- readLines(str_glue("{out_dir}/model.txt"), 4)
        nsample <- str_subset(lines, "Samples") |> str_split_i(" ", 3) |> as.numeric()
        burnin <- str_subset(lines, "Burnin")  |> str_split_i(" ", 3) |> as.numeric()
        burnprop <- burnin / nsample
        thinto <- 1e4L
        priors <- data.table(parameter = names(x),
                             prior = NA_real_)
    }

    x <- list.files(str_glue("{out_dir}"),
                    "trace",
                    full.names = TRUE) |>
        str_subset("combine", negate = TRUE) |>
        str_sort(numeric = TRUE) |>
        map(fread) |>
        map(~ .x[, .SD, .SDcols = !patterns("^L_|Prior|Posterior|Number|log")]) |>
        map(~ setnames(.x, \(x) str_replace(x, "Group effect ", "G_")))

    sildt <- str_chars("sildt")

    walk(x, \(xi) {
        if (FALSE) {
            xi <- copy(x[[1]])
        }

        walk(sildt, \(t1) {
            gt <- str_c("cov_G_", t1, t1)
            et <- str_c("cov_E_", t1, t1)
            pt <- str_c("cov_P_", t1, t1)
            ht <- str_c("h2_", t1, t1)
            if (all(c(gt, et) %in% names(xi))) {
                xi[, (pt) := get(gt) + get(et)]
                xi[, (ht) := get(gt) + get(pt)]
            }
        })

        pnames <- c(str_subset(names(xi), "_P_"),
                    str_subset(names(xi), "h2_"))
        if (length(pnames) > 0) {
            setcolorder(xi, pnames,
                        after = last(str_subset(names(xi), "_[GE]_")),
                        skip_absent = TRUE)
        }
    })

    xs <- x |> map(~ .x[-seq_len(burnin * .N)]) |> rbindlist()
    xs[, state := .I]

    # Write to trace_combine.tsv
    tc <- str_glue("{out_dir}/trace_combine.tsv")
    message(str_glue("- Writing '{tc}'"))
    fwrite(xs, file = tc, sep = "\t")

    xs[, state := NULL]

    xp <- xs[, .(parameter = .SD |> names(),
                 true_val  = NA_real_,
                 mean      = .SD |> map_dbl(mean),
                 median    = .SD |> map_dbl(median),
                 sd        = .SD |> map_dbl(sd),
                 ci95min   = .SD |> map_dbl(quantile, 0.025),
                 ci95max   = .SD |> map_dbl(quantile, 0.975),
                 hdi95min  = .SD |> map(hdi) |> map_dbl("lower"),
                 hdi95max  = .SD |> map(hdi) |> map_dbl("upper"))]

    if (str_detect(dataset, "sim")) {
        tvs <- with(priors, setNames(true_val, parameter))
        xp[parameter %in% names(tvs), true_val := tvs[parameter]]
    }


    # Add in ESS and GR

    # Need to avoid cols with 0 variance
    check <- xs[, map_dbl(.SD, sd)]
    drop_pars <- c("state", names(check[check == 0]))
    cols <- setdiff(xp$parameter, drop_pars)

    # In case we didn't finish or have an uneven number of samples
    min_rows <- map_int(x, nrow) |> min()

    tfs <- x |>
        map(~ .x[seq_len(min_rows), ..cols]) |>
        map(~ mcmc(.x, start = nsample * burnprop, thin = nsample / thinto))

    grd <- as.mcmc.list(tfs) |> gelman.diag() |> _$psrf[, 2]
    ess <- effectiveSize(tfs) |> round()

    xp[parameter %in% cols, `:=`(ESS = ess, GR = grd)]
    xp[, convergence := fcase(ESS >= 500 & GR < 1.05, "***",
                              ESS >= 200 & GR < 1.1, "**",
                              ESS >= 100 & GR < 1.2, "*",
                              default = "")]

    res_file <- str_glue("{res_dir}/{name}.rds")
    if (file.exists(res_file)) {
        tmp <- readRDS(res_file)
        tmp$parameter_estimates <- xp
        saveRDS(tmp, res_file)
    }

    posteriors_tsv <- str_glue("{out_dir}/posterior.tsv")
    message(str_glue("- Writing '{posteriors_tsv}'"))
    fwrite(xp, posteriors_tsv, sep = "\t")
    xp
}

