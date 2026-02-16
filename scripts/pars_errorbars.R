{
    library(data.table)
    library(purrr)
    library(stringr)
    library(ggplot2)
    library(cowplot)

    source("rename_pars.R")
    source("fes_to_vals.R")
}

pars_errorbars <- function(dataset = "fb-test", scens = 0, st_str = "", alt = "") {
    if (FALSE) {
        dataset <- "fb-test"; scens <- 0; st_str = ""; alt <- ""
        dataset <- "sim-base-inf"; scens <- 0; st_str <- "Validating BICI - No. of events"; alt <- "events"
    }

    {
        base_dir <- str_glue("datasets/{dataset}")
        data_dir <- str_glue("{base_dir}/data")
        res_dir  <- str_glue("{base_dir}/results")
        gfx_dir  <- str_glue("{base_dir}/gfx")

        if (!dir.exists(gfx_dir)) {
            message("- mkdir ", gfx_dir)
            dir.create(gfx_dir)
        }
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
        pe[!str_starts(parameter, "Group effect|G_")]
    }) |>
        rbindlist(idcol = "scen", fill = TRUE)

    rf <- str_glue("{res_dir}/{scens_str[[1]]}.rds")
    x[, parameter := rename_bici_pars(parameter)]

    # Extract the widest priors for all parameters
    priors <- map(scens, ~ {
        files <- list.files(res_dir, str_glue("scen-{.x}-"), full.names = TRUE) |>
            str_sort(numeric = TRUE)
        if (is_empty(files)) return(NULL)
        readRDS(files[[1]])$params$priors[, .(parameter, type, val1, val2, true_val)]
    }) |>
        rbindlist(idcol = "scen")

    # Make scenario a factor e.g. c(s2, s4, ...)
    priors[, scen := factor(str_c("s", scens[scen]), levels = str_c("s", scens))]
    # x[, scen := factor(str_c("s", scen), levels = str_c("s", scens))]
    x[, scen := scens_str[scen] |> str_split_i("-", 2) |> str_c("s", x = _) |>
          factor(levels = str_c("s", scens))]

    priors[, `:=`(val1 = min(val1, true_val),
                  val2 = max(val2, true_val)),
           parameter]

    pars <- x[, unique(parameter)] |>
        str_subset("^G_|^Group", negate = TRUE)
    tidy_pars <- setNames(rename_pars(pars), pars)

    x1 <- merge(x[parameter %in% pars],
                priors[parameter %in% pars, .(scen, parameter, type)],
                by = c("scen", "parameter"))

    if ("type" %notin% names(x1)) x1[, type := "uniform"]

    setorder(x1, parameter, median)

    plts <- map(pars, \(par) {
        # i <- 1; par <- pars[[i]]
        y_rng <- priors[parameter == par, c(min(val1), max(val2))]
        y_true <- priors[parameter == par, true_val]
        ymin <- x1[parameter == par, min(hdi95min)]
        ymin <- ymin - 0.1 * abs(ymin)
        ymax <- x1[parameter == par, max(hdi95max)]
        ymax <- ymax + 0.1 * abs(ymax)

        priors2 <- priors[parameter == par, .(scen = as.integer(scen), true_val)]

        ggplot(x1[parameter == par],
                    aes(x = scen, y = median, colour = type)) +
            # geom_boxplot() +
            geom_errorbar(aes(ymin = hdi95min, ymax = hdi95max),
                          position = position_dodge2(),
                          width = 0.5) +
            {if (is_sim)
                geom_segment(data = priors2,
                             aes(x = scen - 0.5, xend = scen + 0.5, y = true_val, yend = true_val),
                             # aes(x = scen, xend = scen, y = true_val, yend = true_val),
                             colour = "green",
                             linewidth = 0.5,
                             linetype = "dashed")} +
            geom_point(position = position_dodge2(width = 0.5),
                       size = 1) +
            geom_hline(yintercept = y_rng[[2]],
                       colour = "grey", linewidth = 0.5, linetype = "dashed") +
            scale_colour_manual(breaks = c("uniform", "inverse", "constant"),
                                values = c("red", "red", "grey40")) +
            scale_x_discrete(drop = FALSE) +
            # scale_y_discrete(limits = ~ range(.x, y_rng)) +
            expand_limits(y = 0) +
            coord_cartesian(ylim = range(0, ymin, ymax)) +
            labs(x = "Scenario",
                 y = "Value",
                 title = tidy_pars[[par]]) +
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
        "infrat", "empty",    "LP_Tr2,Rec", "DP_Tr2,Rec", "RP_Tr2,Rec"
    ) |>
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

    if (str_length(alt) > 0) alt <- str_c("-", alt)
    plt_str <- str_glue("{gfx_dir}/{dataset}-all_hpdi{alt}")

    # ggsave(str_glue("{plt_str}.png"), plt, width = 20, height = 25)
    ggsave(str_glue("{plt_str}.pdf"), plt, width = 20, height = 28.2)
}

if (FALSE) {
    # pars_errorbars("testing", 1, "")
    # pars_errorbars("fb-final", 1:8, "Testing Unlinked vs Linked Traits and FEs")
    # pars_errorbars("fb-final2", 1:4, "Testing Unlinked vs Linked Traits and FEs")
    # pars_errorbars("fb-lp", 1:12, "Testing varying the LP")
    # pars_errorbars("fb-donors", 1:3, "Testing reclassifying Seeder fish")
    # pars_errorbars("fb-simple", 1:6, "Testing if we can fit the LP")
    # pars_errorbars("fb-simple-b", 1:6, "Testing if we can fit the LP with BICI")
    pars_errorbars("sim-test2", 1:5, "Testing coverage")
    pars_errorbars("fb-dp", 1:6, "Testing DPs")
    pars_errorbars("fb-donors", 1:27, "Testing donor reclassification")
    pars_errorbars("fb-test", 0, "Testing BICI on FB data")
    pars_errorbars("fb-test-1e7", 0, "Testing BICI on FB data")
    pars_errorbars("fb-qtest", 0, "Testing BICI on FB data")
    pars_errorbars("sim-base-inf", 0, "Validating BICI")
    pars_errorbars("sim-base-inf", 1:2, "Validating BICI - Base models", "base")
    pars_errorbars("sim-base-inf", 1:12, "Validating BICI - Misspecifying model", "misspecify")
    pars_errorbars("sim-base-inf", c(1:2, 13:20), "Validating BICI - convergence", "conv")
    pars_errorbars("sim-test-inf", 0, "Validating BICI on FB data")
}
