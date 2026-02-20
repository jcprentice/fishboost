{
    library(stringr)
    library(data.table)
    library(purrr)
    library(ggplot2)
    library(ggcorrplot)
}

plot_correlations <- function(dataset = "fb-final", scen = 1, rep = 1) {
    # dataset <- "fb-test"; scen <- 1; rep <- 1;

    message(str_glue("Plotting correlations for '{dataset}/s{scen}-{rep}'"))

    {
        base_dir <- str_glue("datasets/{dataset}")
        data_dir <- str_glue("{base_dir}/data")
        res_dir  <- str_glue("{base_dir}/results")
        gfx_dir  <- str_glue("{base_dir}/gfx")
        cor_dir  <- str_glue("{gfx_dir}/correlations")
    }

    c(gfx_dir, cor_dir) |>
        discard(dir.exists) |>
        walk(~ message("- mkdir ", .x)) |>
        walk(dir.create)

    x <- if (rep == 0) {
        files <- list.files(data_dir, "trace_combine.tsv",
                   recursive = TRUE, full.names = TRUE) |>
            str_sort(numeric = TRUE) |>
            str_subset(str_glue("scen-{scen}-"))
        if (length(files) == 0) {
            return(NULL)
        }
        files |> map(fread) |> rbindlist()
    } else {
        f <- str_glue("{data_dir}/scen-{scen}-{rep}-out/trace_combine.tsv")
        if (!file.exists(f)) {
            return(NULL)
        }
        fread(f)
    }

    params <- if (rep == 0) {
        list.files(res_dir, full.names = TRUE) |>
            str_subset(str_glue("scen-{scen}")) |>
            str_sort(numeric = TRUE) |>
            first() |> readRDS() |> _$params
    } else {
        str_glue("{res_dir}/scen-{scen}-{rep}.rds") |>
            readRDS() |> _$params
    }


    x[, str_subset(names(x), "state|State|Group|^G_") := NULL]

    nx <- str_replace_all(names(x),
                          c("latent_period" = "LP",
                            "detection_period" = "DP",
                            "removal_period" = "RP",
                            "," = "_"))
    setnames(x, nx)

    # Genetic covariances and correlations
    # GVs <- c("cov_G_ss", "cov_G_ii", "cov_G_tt", "r_G_si", "r_G_st", "r_G_it")
    # EVs <- str_replace(GVs, "G", "E")
    GVs <- str_subset(nx, "_G_")
    EVs <- str_subset(nx, "_E_")
    # parameter associated with each fixed effect
    sildt <- c("s", "i", "l", "d", "t")
    fe_par <- setNames(c("beta", "beta", "LP", "DP", "RP"),
                       sildt)
    # FE names
    fes <- c("trial", "donor", "txd", "weight")



    my_theme <- function() {
        theme_classic() +
            theme(axis.title = element_blank(),
                  axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    }

    plts <- map(sildt, \(fe) {
        pars <- c(GVs, fe_par[[fe]], str_c(fes, "_", fe))

        # Some sets don't have weight
        if (!any(str_detect(nx, "weight"))) {
            pars <- str_subset(pars, "^weight", negate = TRUE)
        }

        # Weight may be nested
        if (any(str_detect(nx, "weight1"))) {
            pars <- c(str_subset(pars, "weight", negate = TRUE),
                      str_c("weight", 1:2, "_", fe))
        }

        # This just makes sure that pars actually exist
        pars <- intersect(pars, nx)
        if (is_empty(pars)) return(NULL)

        ggcorrplot(cor(x[, ..pars]), method = "circle") +
            my_theme()

        # datasets/dataset/gfx/dataset-N-cor_l.png
        # ggsave(str_glue("{gfx_dir}/correlations/{dataset}-s{scen}-corr_{fe}.png"),
        #        plt, width = 6, height = 6)
    }) |>
        setNames(sildt)

    description <- params$description |>
        str_split_1(", ") |>
        str_subset("GRM|convergence|coverage", negate = TRUE) |>
        str_flatten_comma() |>
        str_replace_all(c("inf_model 1" = "inf: I = D",
                          "inf_model 2" = "inf: I = 0.1*D",
                          "inf_model 3" = "inf: Don = 0.1*Rec",
                          "inf_model 4" = "inf: Don = p*Rec"))

    title_str <- str_glue("{dataset}/s{scen}-{rep}: {description}")

    plts$all <- ggcorrplot(cor(x),
                           method = "circle",
                           title = title_str) +
        my_theme()

    cov_str <- str_glue("{cor_dir}/{dataset}-s{scen}-{rep}")

    f <- str_glue("{cov_str}-corr_all.pdf")
    message(str_glue("- Plotting '{f}'"))
    ggsave(f, plts$all, width = 12, height = 12)

    GxEpars <- intersect(c(GVs, EVs), nx)
    plts$GxE <- if (length(GxEpars) > 0) {
        ggcorrplot(cor(x[, ..GxEpars]),
                   method = "circle",
                   title = title_str) +
            my_theme()
    }

    f <- str_glue("{cov_str}-corr_GxE.pdf")
    message(str_glue("- Plotting '{f}'"))
    ggsave(f, plts$GxE, width = 6, height = 6)

    pars_G_lp <- str_subset(nx, "cov_G|LP|DP|RP")

    if (length(pars_G_lp) > 0) {
        plts$G_lp <- ggcorrplot(cor(x[, ..pars_G_lp]),
                                method = "circle",
                                title = title_str) +
            my_theme()

        f <- str_glue("{cov_str}-corr_G_lp.pdf")
        message(str_glue("- Plotting '{f}'"))
        ggsave(f, plts$G_lp, width = 12, height = 12)
    }

    plts
}
