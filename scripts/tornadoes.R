{
    library(purrr)
    library(stringr)
    source("plotting/plot_tornadoes.R")
}

tornadoes <- function(dataset = "fb-test", scens = 1:8, combine = TRUE) {

    if (FALSE) {
        dataset <- "sim-test-inf"
        scens <- 1:10
        combine <- TRUE
    }

    gfx_dir <- str_glue("datasets/{dataset}/gfx/tornadoes")
    if (!dir.exists(gfx_dir)) {
        message("- mkdir ", gfx_dir)
        dir.create(gfx_dir, recursive = TRUE)
    }

    plts <- map(scens, \(scen)
                plot_tornadoes(dataset, scen, combine)) |>
        setNames(str_c("s", scens))

    # walk(scens, \(i) {
    #     ggsave(str_glue("{gfx_dir}/{dataset}-s{i}-tornadoes.pdf"),
    #            plts[[i]]$pars, width = 18, height = 12)
    #     ggsave(str_glue("{gfx_dir}/{dataset}-s{i}-tornadoes.png"),
    #            plts[[i]]$pars, width = 18, height = 12)
    # })


    # Create a plot with subplots aligned by trait
    empty <- ggplot() + theme_classic()
    sildt1 <- c("s", "i", "l", "d", "t")
    sildt2 <- str_c(sildt1, sildt1)
    any_non_empty <- function(x) any(x != "empty")

    walk2(plts, names(plts), \(x, scen) {
        # scen <- "s1"; x <- plts[[scen]]
        if (is.null(x)) return()

        x$plots$empty <- empty
        pars <- names(x$plots)

        cov_pars <- c(str_c("cov_G_", sildt2),
                      "r_G_si", "r_G_st", "empty", "empty", "r_G_it",
                      str_c("cov_E_", sildt2),
                      str_c("cov_P_", sildt2))

        model_pars <- c(
            "sigma",  "beta_Tr1", "LP_Tr1,Don", "DP_Tr1,Don", "RP_Tr1,Don",
            "infrat", "empty",    "LP_Tr1,Rec", "DP_Tr1,Rec", "RP_Tr1,Rec",
            "sigma",  "beta_Tr2", "LP_Tr2,Don", "DP_Tr2,Don", "RP_Tr2,Don",
            "infrat", "empty",    "LP_Tr2,Rec", "DP_Tr2,Rec", "RP_Tr2,Rec"
        )

        # Remove repeated sigma and infrat
        beta_in <- str_subset(pars, "beta")
        if (beta_in[[1]] == "beta_Tr2") {
            model_pars[c(1, 6)] <- "empty"
        } else {
            model_pars[c(11, 16)] <- "empty"
        }

        fes <- expand.grid(sildt1,
                           c("trial", "donor", "txd", "weight", "weight1", "weight2")) |>
            rev() |> apply(1, str_flatten, "_")

        plt_names <- c(cov_pars, model_pars, fes) |>
            str_replace_all(c("LP" = "latent_period",
                              "DP" = "detection_period",
                              "RP" = "removal_period"))


        # Some entries like "trial_s" might be missing
        plt_names[plt_names %notin% pars] <- "empty"

        # This clips any rows or columns that are entirely empty
        plt_mat <- matrix(plt_names, nrow = 5)
        plt_mat <- plt_mat[
            which(apply(plt_mat, 1, any_non_empty)),
            which(apply(plt_mat, 2, any_non_empty))
        ]
        plt_names <- as.character(plt_mat)

        pltlst <- with(x$plots, mget(plt_names))

        plt <- plot_grid(x$title_plt,
                         plot_grid(plotlist = pltlst,
                                   ncol = nrow(plt_mat)),
                         ncol = 1,
                         rel_heights = c(0.05, 1))

        plt_str <- str_glue("{gfx_dir}/{dataset}-{scen}-tornadoes")
        ggsave(str_c(plt_str, ".png"), plt, width = 12, height = 12)
        ggsave(str_c(plt_str, ".pdf"), plt, width = 12, height = 12)
        message(str_glue("plotted '{plt_str}'"))
    })
}

# run tests ----
if (FALSE) {
    # tornandoes(dataset = "sim-test2", scens = 1:5, combine = TRUE)
    tornadoes(dataset = "sim-base-inf", scens = 1:10, combine = TRUE)
    tornadoes(dataset = "sim-test-inf", scens = 1:10, combine = TRUE)
}
