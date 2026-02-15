{
    library(data.table)
    library(stringr)
    library(purrr)
    library(coda)
    library(HDInterval)
}

rebuild_sire_posteriors <- function(dataset = "fb-final",
                                    name = "scen-1-1") {

    # dataset <- params$dataset; name <- params$name

    data_dir <- str_glue("datasets/{dataset}/data")
    out_dir  <- str_glue("{data_dir}/{name}-out")
    res_dir  <- str_glue("datasets/{dataset}/results")

    tc_txt <- str_glue("{out_dir}/trace_combine.txt")
    tc     <- str_glue("{out_dir}/trace_combine.tsv")

    if (file.exists(tc_txt)) file.rename(tc_txt, tc)

    x <- fread(tc)
    x[, state := NULL]

    xp <- xs[, .(parameter = .SD |> names(),
                 mean      = .SD |> map_dbl(mean),
                 median    = .SD |> map_dbl(median),
                 ci95min   = .SD |> map_dbl(quantile, 0.025),
                 ci95max   = .SD |> map_dbl(quantile, 0.975),
                 hdi95min  = .SD |> map(hdi) |> map_dbl("lower"),
                 hdi95max  = .SD |> map(hdi) |> map_dbl("upper"))]


    # Add in ESS and GR
    lines <- readLines(str_glue("{out_dir}/model.txt"), 4)

    nsample <- str_subset(lines, "Samples") |> str_split_i(" ", 3) |> as.numeric()
    burnin  <- str_subset(lines, "Burnin")  |> str_split_i(" ", 3) |> as.numeric()

    # Check for txt files and rename to tsv
    txt_files <- list.files(out_dir, "trace.*txt", full.names = TRUE)
    file.rename(txt_files,
                str_replace(txt_files, "txt", "tsv"))

    files <- list.files(out_dir, "trace", full.names = TRUE) |>
        str_subset("extended|combine", negate = TRUE) |>
        str_sort(numeric = TRUE)


    # Need to avoid cols with 0 variance
    check <- x[, map_dbl(.SD, sd)]
    drop_pars <- names(check[check == 0])
    cols <- setdiff(xp$parameter, drop_pars)

    tfs <- map(files, ~ mcmc(fread(.x)[, ..cols]), thin = nsample / 1e4)
    grd <- as.mcmc.list(tfs) |> gelman.diag() |> {\(x) as.data.table(x$psrf)[, 2]}() |> unlist()
    ess <- effectiveSize(tfs) |> round()

    xp[parameter %in% cols, `:=`(ESS = ess, GR = grd)]

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


rebuild_bici_posteriors <- function(dataset = "fb-test",
                                    name = "scen-1-1") {

    # dataset <- params$dataset; name <- params$name

    data_dir <- str_glue("datasets/{dataset}/data")
    res_dir  <- str_glue("datasets/{dataset}/results")
    out_dir  <- str_glue("{data_dir}/{name}-out")

    pfiles <- list.files(str_glue("{out_dir}/output-inf"),
                         "^param_",
                         full.names = TRUE) |>
        str_sort(numeric = TRUE)

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


    change_names <- function(x) {
        x |> str_replace_all(
            c("\\\\" = "", "\\^" = "_",
              "gamma_shape" = "RP_shape",
              "LP" = "latent_period",
              "DP" = "detection_period",
              "RP" = "removal_period",
              "State" = "state",
              "Ω" = "Omega", "ω" = "omega", "μ" = "mu",
              # "^G_" = "Group effect ",
              "cv_G" = "sigma",
              "mu_weight([0-9]?)([sildt])" = "weight\\1_\\2",
              "_gen" = "", "_env" = "",
              "Omega_(.?)g,(.?)g" = "cov_G_\\1\\2",
              "Omega_(.?)e,(.?)e" = "cov_E_\\1\\2",
              "omega_(.?)g,(.?)g" = "r_G_\\1\\2",
              "omega_(.?)e,(.?)e" = "r_E_\\1\\2"))
    }

    x <- pfiles |>
        map(fread) |>
        map(~ .x[, str_subset(names(.SD), "L\\^|N\\^|Prior") := NULL]) |>
        map(~ set_names(.x, change_names(names(.x))))

    xs <- map(x, ~ .x[-seq_len(burnprop * .N)]) |> rbindlist()
    xs[, state := .I]

    # Write to trace_combine.tsv
    tc <- str_glue("{out_dir}/trace_combine.tsv")
    message(str_glue("- Writing '{tc}'"))
    fwrite(xs, file = tc, sep = "\t")
    xs[, state := NULL]

    xp <- xs[, .(parameter = .SD |> names(),
                 mean      = .SD |> map_dbl(mean),
                 median    = .SD |> map_dbl(median),
                 ci95min   = .SD |> map_dbl(quantile, 0.025),
                 ci95max   = .SD |> map_dbl(quantile, 0.975),
                 hdi95min  = .SD |> map(hdi) |> map_dbl("lower"),
                 hdi95max  = .SD |> map(hdi) |> map_dbl("upper"))]

    # Add in ESS and GR

    # Need to avoid cols with 0 variance
    check <- xs[, map_dbl(.SD, sd)]
    drop_pars <- c("state", names(check[check == 0]))
    cols <- setdiff(xp$parameter, drop_pars)

    tfs <- x |>
        map(~ .x[, ..cols]) |>
        map(~ mcmc(.x, start = nsample * burnprop, thin = nsample / thinto))

    grd <- as.mcmc.list(tfs) |> gelman.diag() |> {\(x) as.data.table(x$psrf)[, 2]}() |> unlist()
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
