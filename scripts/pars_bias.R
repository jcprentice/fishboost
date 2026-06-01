{
    library(data.table)
    library(purrr)
    library(stringr)
    library(ggplot2)
    library(cowplot)

    source("rename_pars.R")
    source("fes_to_vals.R")
}

pars_bias <- function(dataset = "fb-test", scens = 0, st_str = "", alt = "", as_grid = TRUE) {
    if (FALSE) {
        dataset <- "fb-test"; scens <- 0; st_str = ""; alt <- ""
        dataset <- "sim-test-inf1"; scens <- 0; st_str <- "Validating BICI"; alt <- ""
        as_grid <- TRUE
    }

    {
        base_dir <- str_glue("datasets/{dataset}")
        res_dir  <- str_glue("{base_dir}/results")
        gfx_dir  <- str_glue("{base_dir}/gfx")

        c(res_dir, gfx_dir) |>
            discard(dir.exists) |>
            walk(~ message(" - mkdir ", .x)) |>
            walk(dir.create)
    }


    is_sim <- str_detect(dataset, "sim")

    scens_str <- list.files(res_dir) |>
        str_remove("\\.rds") |>
        str_sort(numeric = TRUE)

    if (any(scens == 0)) {
        scens <- scens_str |> str_split_i("-", 2) |> unique() |> as.integer()
    } else {
        scens_str <- scens_str |>
            keep(~ .x |> str_split_i("-", 2) |> as.integer() |> is.element(scens))
    }

    # Put scens_str back in the right order
    tmp <- data.table(str = scens_str)
    tmp[, scen := str_split_i(str, "-", 2) |> as.integer(), .I]
    tmp[, pos := match(scen, scens)]
    scens_str <- tmp[order(pos), str]
    rm(tmp)

    x <- map(scens_str, ~ {
        rf <- str_glue("{res_dir}/{.x}.rds")
        pe <- readRDS(rf)$parameter_estimates
        pe[!str_starts(parameter, "G_|Group"),
           .(parameter,
             bias1 = mean - true_val,
             bias2 = (median - true_val) / sd,
             bias3 = median / true_val - 1,
             convergence)]
    }) |>
        rbindlist(idcol = "scen", fill = TRUE)

    x[, scen := scens_str[scen] |> str_split_i("-", 2) |> str_c("s", x = _) |>
          factor(levels = str_c("s", scens))]

    x[, parameter := rename_bici_pars(parameter)]

    pars <- x[, unique(parameter)]
    pretty_pars <- setNames(pretty_names(pars), pars)

    setorder(x, parameter, bias2)

    plts <- map(pars, \(par) {
        # i <- 1; par <- pars[[i]]
        x1 <- x[parameter == par, .(scen, bias = bias2, convergence)]
        mu_x1 <- x1[, .(mu = mean(bias)), scen]
        mu_x1[, scen := as.integer(scen)]

        ggplot(x1, aes(x = scen, y = bias, colour = convergence)) +
            geom_boxplot(fill = "tomato",
                         colour = "black",
                         staplewidth = 0.5,
                         width = 0.3) +
            geom_segment(data = mu_x1,
                         aes(x = scen - 0.4, xend = scen + 0.4,
                             y = mu, yend = mu),
                      colour = "blue",
                      linewidth = 0.5) +
            # geom_point() +
            geom_hline(yintercept = 0, linetype = "dashed") +
            scale_colour_manual(breaks = c("", "*", "**", "***"),
                                values = c("blue3", "green4", "yellow3", "red2")) +
            scale_x_discrete(drop = FALSE) +
            # scale_y_discrete(limits = ~ range(.x, 0, max_bias)) +
            expand_limits(y = 0) +
            labs(x = "Scenario",
                 y = "Bias",
                 title = pretty_pars[[par]]) +
            theme_classic() +
            theme(legend.position = "none")
    }) |> setNames(pars)

    title_plt <- ggplot() +
        labs(title = str_glue("Dataset: '{dataset}'"),
             subtitle = st_str) +
        theme_classic() +
        theme(plot.title = element_text(size = 22),
              plot.subtitle = element_text(size = 16))

    plts$empty <- ggplot() + theme_classic()

    if (as_grid) {

        sildt1 <- c("s", "i", "l", "d", "t")
        sildt2 <- str_c(sildt1, sildt1)
        any_non_empty <- function(x) any(x != "empty")

        cov_pars <- c(str_c("cov_G_", sildt2),
                      "r_G_si", "r_G_st", "empty", "empty", "r_G_it",
                      str_c("cov_E_", sildt2),
                      str_c("cov_P_", sildt2))

        model_pars <- c(
            "sigma",  "beta_Tr1", "LP_Tr1,Don", "DP_Tr1,Don", "RP_Tr1,Don",
            "infrat", "empty",    "LP_Tr1,Rec", "DP_Tr1,Rec", "RP_Tr1,Rec",
            "sigma",  "beta_Tr2", "LP_Tr2,Don", "DP_Tr2,Don", "RP_Tr2,Don",
            "infrat", "empty",    "LP_Tr2,Rec", "DP_Tr2,Rec", "RP_Tr2,Rec") |>
            str_replace_all(c("LP" = "latent_period",
                              "DP" = "detection_period",
                              "RP" = "removal_period"))

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

        plt_names <- c(cov_pars, model_pars, fes)

        # Some entries like "trial_s" might be missing
        plt_names[plt_names %notin% pars] <- "empty"

        # This clips any rows or columns that are entirely empty
        plt_mat <- matrix(plt_names, nrow = 5)
        plt_mat <- plt_mat[
            which(apply(plt_mat, 1, any_non_empty)),
            which(apply(plt_mat, 2, any_non_empty))
        ]
        plt_names <- c(plt_mat)

        pltlst <- with(plts, mget(plt_names))

        plt <- plot_grid(title_plt,
                         plot_grid(plotlist = pltlst,
                                   ncol = nrow(plt_mat),
                                   align = "v"),
                         ncol = 1, rel_heights = c(0.06, 1))
    } else {
        plt <- plot_grid(title_plt,
                         plot_grid(plotlist = plts,
                                   align = "v"),
                         ncol = 1, rel_heights = c(0.06, 1))
    }

    if (str_length(alt) > 0) alt <- str_c("-", alt)
    plt_str <- str_glue("{gfx_dir}/{dataset}-all_bias{alt}")

    message(str_glue("plotted '{plt_str}'"))
    ggsave(str_glue("{plt_str}.png"), plt, width = 20, height = 25)
    ggsave(str_glue("{plt_str}.pdf"), plt, width = 20, height = 25)

    plt
}

if (FALSE) {
    # pars_bias("testing", 1, "")
    # pars_bias("fb-final", 1:8, "Testing Unlinked vs Linked Traits and FEs")
    # pars_bias("fb-final2", 1:4, "Testing Unlinked vs Linked Traits and FEs")
    # pars_bias("fb-lp", 1:12, "Testing varying the LP")
    # pars_bias("fb-donors", 1:3, "Testing reclassifying Seeder fish")
    # pars_bias("fb-simple", 1:6, "Testing if we can fit the LP")
    # pars_bias("fb-simple-b", 1:6, "Testing if we can fit the LP with BICI")
    pars_bias("sim-test2", 1:5, "Testing coverage")
    pars_bias("fb-dp", 1:6, "Testing DPs")
    pars_bias("fb-donors", 1:27, "Testing donor reclassification")
    pars_bias("fb-test", 0, "Testing BICI on FB data")
    pars_bias("fb-test-1e7", 0, "Testing BICI on FB data")
    pars_bias("fb-qtest", 0, "Testing BICI on FB data")
    pars_bias("sim-base-inf", 0, "Validating BICI")
    pars_bias("sim-base-inf", 1:2, "Validating BICI - Base models", "base")
    pars_bias("sim-base-inf", 1:12, "Validating BICI - Misspecifying model", "misspecify")
    # pars_bias("sim-base-inf", c(1:2, 13:20), "Validating BICI - convergence", "conv")
    pars_bias("sim-test-inf", 0, "Validating BICI on Simulated data")
    pars_bias("sim-test-inf", 1:10, "Validating BICI on Simulated data", "tr1")
    pars_bias("sim-test-inf", 11:20, "Validating BICI on Simulated data", "tr12")
}
