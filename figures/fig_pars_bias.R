{
    library(data.table)
    library(stringr)
    library(purrr)
    library(HDInterval)
    library(ggplot2)
    library(cowplot)
    library(ggtext)

    source("rename_pars.R")
    source("figures/theme_natcom.R")
}

fig_pars_bias <- function(dataset = "sim-test-inf1", scens = 0) {
    if (FALSE) {
        dataset <- "sim-test-inf1"; scens <- 0
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
        setup <- res$params$setup |> str_extract("\\d+") |> as.integer()
        pe[str_starts(parameter, "weight_"),
           parameter := str_replace(parameter, "weight_", str_c("weight", setup, "_"))]
        pe[!str_starts(parameter, "G_|Group"),
           .(parameter,
             bias1 = median - true_val,
             bias2 = (median - true_val) / sd)]
    }) |>
        rbindlist(idcol = "scen", fill = TRUE) |>
        _[, scen := scens_str[scen] |> str_split_i("-", 2) |> as.integer()]

    x[, parameter := rename_bici_pars(parameter)]

    # Extract parameters
    pars <- x[, unique(parameter)]
    html_pars <- html_names(pars) |>
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

    l2p <- 1 / ggplot2::.pt

    plts <- map(pars, \(par) {
        if (FALSE) {
            i <- 1
            par <- pars[[i]]
        }
        x1 <- x[parameter == par, .(scen, bias = bias2)]
        x2 <- x1[, .(mu  = mean(bias),
                     min = hdi(bias)[["lower"]],
                     max = hdi(bias)[["upper"]]),
                 scen]
        mu_x1 <- x1[, .(mu = mean(bias)), scen]
        mu_x1[, scen := as.integer(scen)]

        ggplot(x1) +
            geom_errorbar(aes(x = scen, ymin = min, ymax = max),
                          x2,
                          colour = "red",
                          width = 0.2) +
            geom_point(aes(x = scen, y = mu),
                       x2,
                       colour = "red",
                       size = 2 * l2p) +
            # geom_boxplot(aes(x = scen, y = bias),
            #              x1,
            #              fill = "tomato",
            #              colour = "black",
            #              staplewidth = 0.5,
            #              linewidth = 0.35,
            #              width = 0.3,
            #              outliers = FALSE) +
            # geom_segment(aes(x = scen - 0.4, xend = scen + 0.4,
            #                  y = mu, yend = mu),
            #              mu_x1,
            #              colour = "blue",
            #              linewidth = 0.35) +
            geom_hline(yintercept = 0,
                       linetype = "dashed",
                       linewidth = 0.5 * l2p) +
            scale_x_discrete(drop = FALSE) +
            scale_y_continuous(limits = ~ range(.x, 0)) +
            labs(x = NULL,
                 y = NULL,
                 title = html_pars[[par]]) +
            theme_natcom()
    }) |> setNames(pars)

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

    ht <- if (length(trials) == 1) 17 * 2 / 3 else 17

    plt_str <- str_glue("gfx/{dataset}-pars-bias")
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
