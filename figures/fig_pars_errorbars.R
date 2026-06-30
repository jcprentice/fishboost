{
    library(data.table)
    library(stringr)
    library(purrr)
    library(ggplot2)
    library(cowplot)
    library(ggtext)

    source("rename_pars.R")
    source("figures/theme_natcom.R")
}

fig_pars_errorbars <- function(dataset = "fb-test", scens = 0) {
    if (FALSE) {
        dataset <- "fb-test"; scens <- c(1, 2, 7, 8)
        dataset <- "sim-test-inf1"; scens <- 0
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

    x <- map(scens_str, ~ {
        rf <- str_glue("{res_dir}/{.x}.rds")
        res <- readRDS(rf)
        pe <- res$parameter_estimates
        setup <- res$params$setup |> str_extract("\\d+") |> as.integer() |> pluck(1)
        pe[str_starts(parameter, "weight_"),
           parameter := str_replace(parameter, "weight_", str_c("weight", setup, "_"))]
        pe[!str_starts(parameter, "G_|Group"),
           .(parameter, true_val, median, hdi95min, hdi95max, convergence)]
    }) |>
        rbindlist(idcol = "scen", fill = TRUE) |>
        _[, scen := scens_str[scen] |> str_split_i("-", 2) |> as.integer()]

    x[, parameter := rename_bici_pars(parameter)]

    # Extract parameters
    pars <- x[, unique(parameter)]
    html_pars <- html_names(pars) |>
        setNames(pars)

    # Extract the widest priors for all parameters
    priors <- map(scens, \(scen) {
        if (FALSE) {
            scen <- 1
        }
        files <- list.files(res_dir, str_glue("scen-{scen}-"), full.names = TRUE) |>
            str_sort(numeric = TRUE)
        if (is_empty(files)) return(NULL)
        tmp <- readRDS(files[[1]])$params$priors[, .(parameter, type, val1, val2, true_val)]
        tmp[parameter %in% pars]
    }) |>
        rbindlist(idcol = "scen") |>
        _[, scen := scens[scen]]

    # Make scenario a factor
    priors[, scen := factor(scen)]
    x[, scen := factor(scen)]

    # Remove priors for non-existing parameters
    priors <- merge(priors,
                    unique(x[, .(scen, parameter)]),
                    by = c("scen", "parameter"),
                    all.y = TRUE)

    priors[, `:=`(val1 = min(val1, true_val),
                  val2 = max(val2, true_val)),
           parameter]

    x1 <- merge(x, priors[, .(scen, parameter, type)],
                by = c("scen", "parameter"))

    if ("type" %notin% names(x1)) {
        x1[, type := "uniform"]
    }

    setorder(x1, parameter, scen, median)

    # Set lvl names
    lvls <- if ((dataset == "fb-test" && identical(scens, c(1L, 2L, 7L, 8L))) ||
                dataset == "sim-test-inf1") {
        c("SS\nTr1", "MS\nTr1", "SS\nTr1+2", "MS\nTr1+2")
    } else {
        str_c("s", scens)
    }
    setattr(x1$scen, "levels" , lvls)
    setattr(priors$scen, "levels" , lvls)

    plts <- map(pars, \(par) {
        if (FALSE) {
            i <- 1
            par <- pars[[i]]
        }
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

        true_line <- if (is_sim) {
            geom_segment(aes(x = scen - 0.5,xend = scen + 0.5,
                             y = true_val, yend = true_val),
                         priors2,
                         colour = "green",
                         linewidth = 0.5,
                         linetype = "dashed")
        }


        l2p <- 1 / ggplot2::.pt

        ggplot(x1[parameter == par],
                    aes(x = scen, y = median,
                        group = scen, colour = convergence)) +
            geom_errorbar(aes(ymin = hdi95min, ymax = hdi95max),
                          position = position_dodge2(),
                          width = 0.5,
                          linewidth = 1 * l2p) +
            geom_point(position = position_dodge2(width = 0.5),
                       size = 0.5) +
            true_line +
            scale_colour_manual(breaks = conv_breaks,
                                values = conv_values) +
            scale_x_discrete(drop = FALSE) +
            scale_y_continuous(limits = ~ range(0, ymin, ymax)) +
            labs(x = NULL,
                 y = NULL,
                 title = html_pars[[par]]) +
            theme_natcom()
                  # axis.text.x = element_text(angle = 0, hjust = 1))
    }) |> setNames(pars)


    # May need to rename weights if only 1 trial
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

    plt <- plot_grid(plotlist = plts[plt_names],
                     ncol = 6, align = "v")
    plt

    ht <- if (length(trials) == 1) 12 else 17

    plot_str <- str_glue("gfx/{dataset}-pars-errorbars")
    ggsave(str_glue("{plot_str}.pdf"), plt,
           width = 18.3, height = ht, units = "cm")
    ggsave(str_glue("{plot_str}.png"), plt,
           width = 18.3, height = ht, units = "cm", dpi = "print")

    plt
}

if (FALSE) {
    plt <- fig_pars_errorbars("fb-test", c(1, 2, 7, 8))
    plt <- fig_pars_errorbars("sim-test-inf2", 0)
    plt <- fig_pars_errorbars("sim-test-inf2", 0)
    plt
}
