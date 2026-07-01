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

pars_bias <- function(dataset = "fb-test",
                      scens = 0,
                      st_str = "",
                      alt = "",
                      plt_shape = "traits",
                      output = "pdf") {
    if (FALSE) {
        dataset <- "fb-test"; scens <- 0; st_str = ""; alt <- ""
        dataset <- "sim-test-inf1"; scens <- 0; st_str <- "Validating BICI"; alt <- ""
        plt_shape <- "traits"; output <- "pdf"
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

    l2p <- 1 / ggplot2::.pt

    plts <- map(pars, \(par) {
        if (FALSE) {
            i <- 1; par <- pars[[i]]
        }
        x1 <- x[parameter == par, .(scen, bias = bias2)]
        x2 <- x1[, .(mu  = mean(bias),
                     min = hdi(bias)[["lower"]],
                     max = hdi(bias)[["upper"]]),
                 scen]
        mu_x1 <- x1[, .(mu = mean(bias)), scen]
        mu_x1[, scen := as.integer(scen)]

        ggplot() +
            geom_errorbar(aes(x = scen, ymin = min, ymax = max),
                          x2,
                          colour = "red",
                          width = 0.2) +
            geom_point(aes(x = scen, y = mu),
                       x2,
                       colour = "red",
                       size = 2 * l2p) +
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

    pmat <- get_plot_matrix(pars, plt_shape)

    plt <- plot_grid(title_plt,
                     plot_grid(plotlist = plts[pmat$plt_names],
                               nrow = pmat$nr,
                               ncol = pmat$nc,
                               align = "hv"),
                     ncol = 1,
                     rel_heights = c(0.75, pmat$nr))

    if (str_length(alt) > 0) alt <- str_c("-", alt)

    plt_str <- str_glue("{gfx_dir}/{dataset}-all_bias{alt}")

    message(str_glue("plotted '{plt_str}'"))

    walk(output, \(op) {
        ggsave(str_glue("{plt_str}.{op}"), plt,
               width = 4 * pmat$nc,
               height = 3 * (pmat$nr + 0.75))
    })

    plt
}

if (FALSE) {
    pars_bias("sim-test-inf1", 0, "Validating BICI on Simulated data")
    pars_bias("sim-test-inf2", 0, "Validating BICI on Simulated data")
    pars_bias("sim-events", 0, "Testing BICI's handling of events")
}
