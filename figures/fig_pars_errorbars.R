{
    library(data.table)
    library(stringr)
    library(purrr)
    library(ggplot2)
    library(cowplot)
    library(ggtext)

    source("rename_pars.R")
}

fig_pars_errorbars <- function() {
    dataset <- "fb-test"
    scens <- c(1, 2, 7, 8)

    {
        base_dir <- str_glue("datasets/{dataset}")
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
        scens <- scens_str |> str_split_i("-", 2) |> unique() |> as.integer()
    } else {
        scens_str <- scens_str |>
            keep(~ .x |> str_split_i("-", 2) |> as.integer() |> is.element(scens))
    }

    # Put scens_str back into the order specified by scens
    tmp <- data.table(str = scens_str)
    tmp[, scen := str_split_i(str, "-", 2) |> as.integer(), .I]
    tmp[, pos := match(scen, scens)]
    scens_str <- tmp[order(pos), str]
    rm(tmp)

    x <- map(scens_str, ~ {
        rf <- str_glue("{res_dir}/{.x}.rds")
        res <- readRDS(rf)
        pe <- res$parameter_estimates
        setup <- res$params$setup |> str_extract("\\d+") |> as.integer()
        pe[str_starts(parameter, "weight_"),
           parameter := str_replace(parameter, "weight_", str_c("weight", setup, "_"))]
        pe[!str_starts(parameter, "G_|Group"),
                       .(parameter, true_val, median, hdi95min, hdi95max, convergence)]
    }) |>
        rbindlist(idcol = "scen", fill = TRUE) |>
        _[, scen := scens[scen]]

    rf <- str_glue("{res_dir}/{scens_str[[1]]}.rds")
    x[, parameter := rename_bici_pars(parameter)]

    # Extract the widest priors for all parameters
    priors <- map(scens, ~ {
        files <- list.files(res_dir, str_glue("scen-{.x}-"), full.names = TRUE) |>
            str_sort(numeric = TRUE)
        if (is_empty(files)) return(NULL)
        readRDS(files[[1]])$params$priors[, .(parameter, type, val1, val2, true_val)]
    }) |>
        rbindlist(idcol = "scen") |>
        _[, scen := scens[scen]]

    # Make scenario a factor
    priors[, scen := factor(scen)]
    x[, scen := factor(scen)]

    priors[, `:=`(val1 = min(val1, true_val),
                  val2 = max(val2, true_val)),
           parameter]

    pars <- x[, unique(parameter)]
    html_pars <- html_names(pars) |>
        # str_replace_all(c(" Gen" = "<sub>A</sub>",
        #                   " Env" = "<sub>E</sub>",
        #                   "Period " = "Period<br>",
        #                   "Trial (.) \\(" = "(Trial \\1, ")) |>
        setNames(pars)

    x1 <- merge(x[parameter %in% pars],
                priors[parameter %in% pars, .(scen, parameter, type)],
                by = c("scen", "parameter"))

    if ("type" %notin% names(x1)) {
        x1[, type := "uniform"]
    }

    setorder(x1, parameter, scen, median)

    # Set lvl names
    lvls <- expand.grid(c("SS", "MS"),
                    c("Tr1", "Tr1+2")) |>
            apply(1, str_flatten, collapse = "\n")
    setattr(x1$scen, "levels" , lvls)
    setattr(priors$scen, "levels" , lvls)

    plts <- map(pars, \(par) {
        # i <- 1; par <- pars[[i]]
        y_rng <- priors[parameter == par, c(min(val1), max(val2))]
        y_true <- priors[parameter == par, true_val]
        ymin <- x1[parameter == par, min(hdi95min)]
        ymin <- ymin - 0.1 * abs(ymin)
        ymax <- x1[parameter == par, max(hdi95max)]
        ymax <- ymax + 0.1 * abs(ymax)

        priors2 <- priors[parameter == par, .(scen = as.integer(scen), true_val)]

        # Colour = type vs colour = convergence

        conv_breaks <- c("", "*", "**", "***")
        # conv_cols <- c("blue3", "green4", "yellow3","red2")
        conv_values <- c("blue3", "blue3", "red2","red2")

        ggplot(x1[parameter == par],
                    aes(x = scen, y = median,
                        group = scen, colour = convergence)) +
            # geom_boxplot() +
            geom_errorbar(aes(ymin = hdi95min, ymax = hdi95max),
                          position = position_dodge2(),
                          width = 0.5,
                          linewidth = 0.35) +
            geom_point(position = position_dodge2(width = 0.5),
                       size = 0.5) +
            # geom_hline(yintercept = y_rng[[2]],
            #            colour = "grey", linewidth = 0.5, linetype = "dashed") +
            # scale_colour_manual(breaks = c("uniform", "inverse", "constant"),
            #                     values = c("red", "red", "grey40")) +
            scale_colour_manual(breaks = conv_breaks,
                                values = conv_values) +
            scale_x_discrete(drop = FALSE) +
            # scale_y_discrete(limits = ~ range(.x, y_rng)) +
            expand_limits(y = 0) +
            coord_cartesian(ylim = range(0, ymin, ymax)) +
            labs(x = NULL,
                 y = NULL,
                 title = html_pars[[par]]) +
            theme_classic() +
            theme(legend.position = "none",
                  text = element_text(size = 5),
                  plot.title = element_markdown(size = 7),
                  axis.line = element_line(linewidth = 0.35),
                  axis.ticks = element_line(linewidth = 0.2))
                  # axis.text.x = element_text(angle = 0, hjust = 1))
    }) |> setNames(pars)

    pltnames <- c("cov_G_ss", "cov_G_ii", "cov_G_tt", "r_G_si", "r_G_st", "r_G_it",
                  "cov_E_ss", "cov_E_ii", "cov_E_tt", "r_E_si", "r_E_st", "r_E_it",
                  "LP_Tr1,Don", "DP_Tr1,Don", "RP_Tr1,Don",
                  "weight1_s", "weight1_i", "weight1_t",
                  "LP_Tr1,Rec", "DP_Tr1,Rec", "RP_Tr1,Rec",
                  "weight2_s", "weight2_i", "weight2_t",
                  "LP_Tr2,Don", "DP_Tr2,Don", "RP_Tr2,Don",
                  "beta_Tr1", "beta_Tr2", "infrat",
                  "LP_Tr2,Rec", "DP_Tr2,Rec", "RP_Tr2,Rec",
                  "sigma") |>
        str_replace_all(c("LP" = "latent_period",
                          "DP" = "detection_period",
                          "RP" = "removal_period"))


    plt <- plot_grid(plotlist = plts[pltnames],
                     ncol = 6, align = "v")
    plt

    ggsave("gfx/fb-pars-errorbars.png", plt,
           width = 18.3, height = 17, units = "cm")
    ggsave("gfx/fb-pars-errorbars.pdf", plt,
           width = 18.3, height = 17, units = "cm", dpi = "print")

    plt
}

fig_pars_errorbars()
