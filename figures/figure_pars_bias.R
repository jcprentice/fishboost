{
    library(data.table)
    library(stringr)
    library(purrr)
    library(ggplot2)
    library(cowplot)
    library(ggtext)

    source("rename_pars.R")
    source("fes_to_vals.R")
}

figure_pars_bias <- function(dataset = "sim-test-inf2", scens = 0) {
    if (FALSE) {
        dataset <- "sim-test-inf2"; scens <- 0
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
        # rf <- str_glue("{res_dir}/{scens_str[[1]]}.rds")
        res <- readRDS(rf)
        pe <- res$parameter_estimates
        setup <- res$params$setup |> str_extract("\\d+") |> as.integer()
        pe[str_starts(parameter, "weight_"),
           parameter := str_replace(parameter, "weight_", str_c("weight", setup, "_"))]
        pe[!str_starts(parameter, "G_|Group"),
           .(parameter,
             desc = res$params$description,
             bias1 = mean - true_val,
             bias2 = (median - true_val) / sd,
             bias3 = median / true_val - 1,
             convergence)]
    }) |>
        rbindlist(idcol = "scen", fill = TRUE) |>
        _[, scen := scens_str[scen] |> str_split_i("-", 2) |> as.integer()]

    descriptions <- x[, unique(desc)]
    x[, desc := NULL]

    x[, parameter := rename_bici_pars(parameter)]

    pars <- x[, unique(parameter)]
    html_pars <- html_names(pars) |>
        # str_replace_all(c(" Gen" = "<sub>A</sub>",
        #                   " Env" = "<sub>E</sub>",
        #                   "Period " = "Period<br>",
        #                   "Trial (.) \\(" = "(Trial \\1, ")) |>
        setNames(pars)

    setorder(x, parameter, bias2)

    # Set lvl names
    lvls <- if (dataset == "sim-test-inf1") {
        c("SS\nTr1", "MS\nTr1", "SS\nTr1+2", "MS\nTr1+2")
    } else if (dataset == "sim-test-inf2") {
        str_c("s", scens)
    } else {
        map_chr(descriptions, ~ {
            str_c(if (str_detect(.x, "GEV SIT,")) "SS" else "MS",
                  "\n",
                  fcase(str_detect(.x, "FB_1_"), "Tr1",
                        str_detect(.x, "FB_2_"), "Tr2",
                        default = "Tr1+2"))
        })
    }
    x[, scen := factor(scen)]
    setattr(x$scen, "levels" , lvls)

    plts <- map(pars, \(par) {
        if (FALSE) {
            i <- 1
            par <- pars[[i]]
        }
        x1 <- x[parameter == par, .(scen, bias = bias2, convergence)]
        mu_x1 <- x1[, .(mu = mean(bias)), scen]
        mu_x1[, scen := as.integer(scen)]

        ggplot(x1) +
            geom_boxplot(aes(x = scen, y = bias),
                         fill = "tomato",
                         colour = "black",
                         staplewidth = 0.5,
                         linewidth = 0.35,
                         width = 0.3,
                         outliers = FALSE) +
            geom_segment(data = mu_x1,
                         aes(x = scen - 0.4, xend = scen + 0.4,
                             y = mu, yend = mu),
                      colour = "blue",
                      linewidth = 0.35) +
            # geom_point() +
            geom_hline(yintercept = 0,
                       linetype = "dashed",
                       linewidth = 0.35) +
            scale_x_discrete(drop = FALSE) +
            # scale_y_discrete(limits = ~ range(.x, 0, max_bias)) +
            expand_limits(y = 0) +
            labs(x = NULL,
                 y = NULL,
                 title = html_pars[[par]]) +
            theme_classic() +
            theme(legend.position = "none",
                  text = element_text(size = 5),
                  plot.title = element_markdown(size = 7),
                  axis.line = element_line(linewidth = 0.35),
                  axis.ticks = element_line(linewidth = 0.2))
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
                          "RP" = "removal_period")) |>
        intersect(names(plts))


    plt <- plot_grid(plotlist = plts[pltnames],
                     ncol = 6, align = "v")
    plt

    ggsave("gfx/fb-pars-bias.png", plt,
           width = 18.3, height = 17, units = "cm")
    ggsave("gfx/fb-pars-bias.pdf", plt,
           width = 18.3, height = 17, units = "cm", dpi = "print")

    plt
}

if (FALSE) {
    figure_pars_bias("sim-test-inf1")
    figure_pars_bias("sim-test-inf2")
}
