{
    library(data.table)
    library(stringr)
    library(purrr)
    library(HDInterval)
    library(ggplot2)
    library(cowplot)
    library(ggtext)

    safe_source("rename_pars.R")
    safe_source("figures/theme_natcom.R")
}

fig_pars_zcontract <- function(dataset = "sim-test-inf1", scens = 0) {
    if (FALSE) {
        dataset <- "sim-test-inf1"; scens <- 0
        dataset <- "sim-test-inf2"; scens <- 0
    }

    {
        base_dir <- str_glue("datasets/{dataset}")
        data_dir <- str_glue("{base_dir}/data")
        res_dir  <- str_glue("{base_dir}/results")
        gfx_dir  <- str_glue("{base_dir}/gfx")

        c(res_dir, gfx_dir) |>
            discard(dir.exists) |>
            walk(~ message(" - mkdir ", .x)) |>
            walk(dir.create)
    }

    scens_str <- list.files(res_dir) |>
        str_remove("\\.rds") |>
        str_sort(numeric = TRUE)

    if (any(scens == 0)) {
        scens <- scens_str |> str_split_i("-", 2) |> unique()
    } else {
        scens_str <- scens_str |>
            keep(~ .x |> str_split_i("-", 2) |> as.integer() |> is.element(scens))
    }
    scens <- as.integer(scens)

    # Put scens_str back into the order specified by scens
    tmp <- data.table(str = scens_str)
    tmp[, scen := str_split_i(str, "-", 2) |> as.integer(), .I]
    tmp[, pos := match(scen, scens)]
    scens_str <- tmp[order(pos), str]
    rm(tmp)

    get_scen <- function(ss) str_split_i(ss, "-", 2) |> as.integer()
    get_rep  <- function(ss) str_split_i(ss, "-", 3) |> as.integer()

    posteriors <- map(scens_str, \(ss) {
        if (FALSE) {
            ss <- "scen-1-1"
        }
        f <- str_glue("{data_dir}/{ss}-out/trace_combine.tsv")
        if (!file.exists(f)) {
            return(NULL)
        }
        tc <- fread(f)

        if (any(str_detect(names(tc), "weight_"))) {
            f <- str_glue("{res_dir}/{ss}.rds")
            setup <- readRDS(f)$params$setup |> str_extract("fb_(\\d+)", 1)
            setnames(tc, \(s) str_replace(s, "weight_", str_glue("weight{setup}_")))
        }

        tc[, names(.SD) := NULL, .SDcols = patterns("state|^G_")]
        tc1 <- melt(tc, measure.vars = names(tc),
                    variable.name = "parameter",
                    variable.factor = FALSE)
        tc2 <- tc1[, .(post_mean = mean(value),
                       post_var = var(value)),
                   parameter]

    }) |>
        rbindlist(idcol = "id", fill = TRUE)

    priors <- map(scens_str, \(ss) {
        if (FALSE) {
            ss <- "scen-1-1"
        }
        rf <- str_glue("{res_dir}/{ss}.rds")
        readRDS(rf)$params$priors |>
            _[, .(parameter, true_val, type, val1, val2)]
    }) |>
        rbindlist(idcol = "id", fill = TRUE)

    inv_var <- function(a, b) {
        (b^2 - a^2) / (2 * log(b / a)) -
            ((b - a) / log(b / a))^2
    }

    {
        priors[type == "uniform", prior_var := (val2 - val1)^2 / 12]
        priors[type == "inverse", prior_var := inv_var(val1, val2), .I]
        priors[str_starts(parameter, "cov_"), prior_var := 4 * (1 - 2 / pi)]
        priors[str_starts(parameter, "r_"), prior_var := 1 / (2 * 1.2 + 1)]
        priors[, c("type", "val1", "val2") := NULL]
    }

    x <- merge(posteriors, priors, by = c("id", "parameter"))
    x[, rank := match(parameter, param_order)]
    setorder(x, id, rank)
    x[, rank := NULL]

    # Add scen and rep columns
    x[, `:=`(scen = scens_str[id] |> str_split_i("-", 2) |> as.integer(),
             rep  = scens_str[id] |> str_split_i("-", 3) |> as.integer())]
    setcolorder(x, c("id", "scen", "rep"))

    pars <- intersect(param_order, x$parameter)

    x[, rank := match(parameter, param_order)]
    setorder(x, id, scen, rep, rank)
    x[, rank := NULL]

    x[, `:=`(zscore = (post_mean - true_val) / sqrt(post_var),
             contraction = 1 - post_var / prior_var)]

    # Set lvl names
    lvls <- if (dataset == "sim-test-inf1") {
        c("SSE Trial 1", "MSE Trial 1", "SSE Trials 1+2", "MSE Trials 1+2")
    } else if (dataset == "sim-test-inf2") {
        str_c("s", scens)
    } else {
        map_chr(scens, \(s) {
            rf <- list.files(res_dir, str_glue("scen-{s}-"),
                       full.names = TRUE) |>
                str_sort(numeric = TRUE) |>
                first() |>
                readRDS()
            d <- rf$params$description
            str_c(if (str_detect(d, "(IEs|GEV) SIT,")) "SSE " else "MSE ",
                  fcase(str_detect(d, "FB_1_"), "Trial 1",
                        str_detect(d, "FB_2_"), "Trial 2",
                        default = "Trials 1+2"))
        })
    }
    x[, scen := factor(scen)]
    setattr(x$scen, "levels" , lvls)

    l2p <- 1 / ggplot2::.pt

    plts <- map(pars, \(par) {
        if (FALSE) {
            i <- 1
            par <- pars[[i]]
        }
        x1 <- x[parameter == par, .(scen, zscore, contraction)]

        ggplot(x1) +
            aes(x = contraction,
                y = zscore,
                colour = scen) +
            geom_hline(yintercept = 0,
                       linetype = "dashed") +
            geom_point() +
            scale_x_continuous(limits = ~ range(., 0, 1)) +
            labs(x = "Contraction",
                 y = "z-score",
                 title = html_names(par)) +
            # facet_wrap(. ~ scen) +
            theme_bw() +
            theme(plot.title = element_markdown())
    }) |>
        setNames(pars)

    trials <- str_extract(pars, "beta_Tr(.)", group = 1) |>
        discard(is.na) |> unique() |> as.integer() |> sort()

    plt_names <- c(
        "cov_G_ss", "cov_G_ii", "cov_G_tt", "r_G_si", "r_G_st", "r_G_it",
        "cov_E_ss", "cov_E_ii", "cov_E_tt", "r_E_si", "r_E_st", "r_E_it",
        if (identical(trials, 1L)) {
            c("LP_Tr1,Don", "DP_Tr1,Don", "RP_Tr1,Don", "weight1_s", "weight1_i", "weight1_t",
              "LP_Tr1,Rec", "DP_Tr1,Rec", "RP_Tr1,Rec", "beta_Tr1",  "infrat",    "sigma")
        } else if (identical(trials, 2L)) {
            c("LP_Tr2,Don", "DP_Tr2,Don", "RP_Tr2,Don", "weight2_s", "weight2_i", "weight2_t",
              "LP_Tr2,Rec", "DP_Tr2,Rec", "RP_Tr2,Rec", "beta_Tr2",  "infrat",    "sigma")
        } else {
            c("LP_Tr1,Don", "DP_Tr1,Don", "RP_Tr1,Don", "weight1_s", "weight1_i", "weight1_t",
              "LP_Tr1,Rec", "DP_Tr1,Rec", "RP_Tr1,Rec", "weight2_s", "weight2_i", "weight2_t",
              "LP_Tr2,Don", "DP_Tr2,Don", "RP_Tr2,Don", "beta_Tr1",  "beta_Tr2",  "infrat",
              "LP_Tr2,Rec", "DP_Tr2,Rec", "RP_Tr2,Rec", "sigma")
        }
    )

    legend <- get_legend(
        # create some space to the left of the legend
        plts[[1]] +
            theme(legend.position = "bottom",
                  legend.box.margin = margin(0, 0, 0, 12))
    )

    plts2 <- map(plts[plt_names], \(p) {
        p + theme(legend.position = "none")
    })

    plt <- plot_grid(plotlist = plts2,
                     ncol = 6, align = "v")
    plt

    plt2 <- plot_grid(plt, legend,
                      ncol = 1,
                      rel_heights = c(1, 0.1))
    plt2

    ht <- if (length(trials) == 1) 17 * 2 / 3 else 17

    plt_str <- str_glue("gfx/{dataset}-pars-zscore")
    ggsave(str_glue("{plt_str}.pdf"), plt,
           width = 18.3, height = ht, units = "cm", dpi = "print")
    ggsave(str_glue("{plt_str}.png"), plt,
           width = 18.3, height = ht, units = "cm", dpi = "print")

    plt
}

if (FALSE) {
    fig_pars_bias("sim-test-inf1")
    fig_pars_bias("sim-test-inf2")
}
