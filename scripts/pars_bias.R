{
    library(data.table)
    library(purrr)
    library(stringr)
    library(ggplot2)
    library(ggtext)
    library(cowplot)

    source("utils.R")
    source("rename_pars.R")
    source("fes_to_vals.R")
    source("figures/theme_natcom.R")
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
             bias3 = median / true_val - 1)]
    }) |>
        rbindlist(idcol = "scen", fill = TRUE)

    x[, scen := scens_str[scen] |> str_split_i("-", 2) |> str_c("s", x = _) |>
          factor(levels = str_c("s", scens))]

    x[, parameter := rename_bici_pars(parameter)]

    pars <- x[, unique(parameter)]
    html_pars <- setNames(html_names(pars), pars)

    setorder(x, parameter, bias2)

    plts <- map(pars, \(par) {
        if (FALSE) {
            i <- 1; par <- pars[[i]]
        }
        x1 <- x[parameter == par, .(scen, bias = bias2)]
        x2 <- x1[, .(mu  = mean(bias),
                     min = hdi(bias)[["lower"]],
                     max = hdi(bias)[["upper"]])]
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

    title_plt <- ggplot() +
        labs(title = str_glue("Dataset: '{dataset}'"),
             subtitle = st_str) +
        theme_classic() +
        theme(plot.title = element_text(size = 22),
              plot.subtitle = element_text(size = 16))

    plts$empty <- ggplot() + theme_classic()

    if (as_grid) {

        sildt1 <- str_chars("sildt")
        sildt2 <- str_c(sildt1, sildt1)
        any_non_empty <- function(x) any(x != "empty")

        cov_pars <- c(str_c("cov_G_", sildt2),
                      "r_G_si", "r_G_st", "empty", "empty", "r_G_it",
                      str_c("cov_E_", sildt2))
        # cov_pars <- c(cov_pars, str_c("cov_P_", sildt2), str_c("h2_", sildt2))

        model_pars <- c(
            "sigma",  "beta_Tr1", "LP_Tr1,Don", "DP_Tr1,Don", "RP_Tr1,Don",
            "infrat", "empty",    "LP_Tr1,Rec", "DP_Tr1,Rec", "RP_Tr1,Rec",
            "sigma",  "beta_Tr2", "LP_Tr2,Don", "DP_Tr2,Don", "RP_Tr2,Don",
            "infrat", "empty",    "LP_Tr2,Rec", "DP_Tr2,Rec", "RP_Tr2,Rec")

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
        plt_mat <- matrix(plt_names, nrow = 5, byrow = TRUE)
        plt_mat <- plt_mat[
            which(apply(plt_mat, 1, any_non_empty)),
            which(apply(plt_mat, 2, any_non_empty))
        ]
        plt_names <- c(t(plt_mat))

        pltlst <- with(plts, mget(plt_names))

        nc <- ncol(plt_mat)
        nr <- nrow(plt_mat)

        plt <- plot_grid(title_plt,
                         plot_grid(plotlist = pltlst,
                                   nrow = nr, ncol = nc,
                                   byrow = TRUE, align = "v"),
                         ncol = 1, rel_heights = c(0.75 / nr, 1))
    } else {
        nc <- 5
        nr <- ceiling(length(plts) / nc)

        plt <- plot_grid(title_plt,
                         plot_grid(plotlist = plts,
                                   nrow = nr, ncol = nc,
                                   byrow = TRUE, align = "v"),
                         ncol = 1, rel_heights = c(0.75 / nr, 1))
    }

    if (str_length(alt) > 0) alt <- str_c("-", alt)

    plt_str <- str_glue("{gfx_dir}/{dataset}-all_bias{alt}")

    message(str_glue("plotted '{plt_str}'"))
    ggsave(str_glue("{plt_str}.png"), plt, width = 4 * nc, height = 3 * (nr + 0.75))
    ggsave(str_glue("{plt_str}.pdf"), plt, width = 4 * nc, height = 3 * (nr + 0.75))

    plt
}

if (FALSE) {
    pars_bias("sim-test-inf1", 0, "Validating BICI on Simulated data")
    pars_bias("sim-test-inf2", 0, "Validating BICI on Simulated data")
    pars_bias("sim-events", 0, "Testing BICI's handling of events")
}
