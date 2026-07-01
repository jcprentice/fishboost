{
    library(data.table)
    library(purrr)
    library(stringr)
    library(ggplot2)
    library(ggtext)
    library(cowplot)

    source("utils.R")
    source("get_plot_matrix.R")
    source("rename_pars.R")
    source("fes_to_vals.R")
}

pars_errorbars <- function(dataset = "fb-test",
                           scens = 0,
                           st_str = "",
                           alt = "",
                           plt_shape = "traits",
                           output = "pdf") {
    if (FALSE) {
        dataset <- "fb-test"; scens <- 0; st_str = ""
        dataset <- "sim-test-inf1"; scens <- 0; st_str <- "Validating BICI"
        alt <- ""; plt_shape <- "traits"; output <- "pdf"
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

    # Put scens_str back into the order specified by scens
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

    # Extract parameters
    pars <- x[, unique(parameter)]
    html_pars <- setNames(html_names(pars), pars)


    # Extract the widest priors for all parameters
    priors <- map(scens, ~ {
        files <- list.files(res_dir, str_glue("scen-{.x}-"), full.names = TRUE) |>
            str_sort(numeric = TRUE)
        if (is_empty(files)) return(NULL)
        tmp <- readRDS(files[[1]])$params$priors[, .(parameter, type, val1, val2, true_val)]
        tmp[parameter %in% pars]
    }) |>
        rbindlist(idcol = "scen")

    # Make scenario a factor e.g. c(s2, s4, ...)
    priors[, scen := factor(str_c("s", scens[scen]), levels = str_c("s", scens))]
    # x[, scen := factor(str_c("s", scen), levels = str_c("s", scens))]
    x[, scen := scens_str[scen] |> str_split_i("-", 2) |> str_c("s", x = _) |>
          factor(levels = str_c("s", scens))]

    # Remove priors for non-existing parameters
    priors <- merge(priors,
                    unique(x[, .(scen, parameter)]),
                    by = c("scen", "parameter"),
                    all.y = TRUE)

    priors[, `:=`(val1 = min(val1, true_val),
                  val2 = max(val2, true_val)),
           parameter]
    priors[str_starts(parameter, "h2"), `:=`(val1 = 0, val2 = 1)]

    x1 <- merge(x, priors[, .(scen, parameter, type)],
                by = c("scen", "parameter"))

    if ("type" %notin% names(x1)) {
        x1[, type := "uniform"]
    }

    setorder(x1, parameter, scen, median)

    # Set lvl names
    lvls <- if (dataset == "fb-test") {
        c("SS End\nTr1", "MS End\nTr1", "No Var\nTr1",
          "SS End\nTr2", "MS End\nTr2", "No Var\nTr2",
          "SS End\nTr1+2", "MS End\nTr1+2", "No Var\nTr1+2")
    } else if (dataset == "sim-test-inf1") {
        c("SS End\nTr1", "MS End\nTr1", "SS End\nTr1+2", "MS End\nTr1+2")
    } else if (dataset == "sim-test-inf2") {
        c("Correct model",
          "Overfitting MSE to No Vars",
          "Overfitting MSE to ST",
          "Testing h2=0",
          "Underfitting MSE to Cors=0",
          "Underfitting MSE to All vars",
          "Underfitting none to SSE")
    } else {
        levels(x1$scen)
    }
    setattr(x1$scen, "levels" , lvls)
    setattr(priors$scen, "levels" , lvls)


    # Sort out phenotypic variance and heritability
    x1[str_detect(parameter, "_P_|h2_"), convergence := "NA"]
    walk(c("ss", "ii", "tt"), \(xx) {
       G <- priors[str_detect(parameter, str_c("_G_", xx)), true_val]
       E <- priors[str_detect(parameter, str_c("_E_", xx)), true_val]
       P <- G + E
       h2 <- G / P
       priors[str_detect(parameter, str_c("_P_", xx)), true_val := P]
       priors[str_detect(parameter, str_c("h2_", xx)), true_val := h2]
    })



    plts <- map(pars, \(par) {
        if (FALSE) {
            i <- 1; par <- pars[[i]]
        }
        y_rng <- priors[parameter == par, c(min(val1), max(val2))]
        y_true <- priors[parameter == par, true_val]
        ymin <- x1[parameter == par, min(hdi95min)]
        ymin <- ymin - 0.1 * abs(ymin)
        ymax <- x1[parameter == par, max(hdi95max)]
        ymax <- ymax + 0.1 * abs(ymax)

        priors2 <- priors[parameter == par, .(scen = as.integer(scen), true_val)]

        # Colour = type vs colour = convergence

        conv_breaks <- c("", "*", "**", "***", "NA")
        # conv_cols <- c("blue3", "green4", "yellow3","red2")
        conv_values <- c("blue3", "blue3", "red2","red2", "grey30")

        true_line <- if (is_sim) {
            geom_segment(aes(x = scen - 0.5, xend = scen + 0.5,
                             y = true_val, yend = true_val),
                         priors2,
                         colour = "green",
                         linewidth = 0.5 / .pt,
                         linetype = "dashed")
        }

        ggplot(x1[parameter == par]) +
            aes(x = scen, y = median, group = scen, colour = convergence) +
            # geom_boxplot() +
            geom_errorbar(aes(ymin = hdi95min, ymax = hdi95max),
                          position = position_dodge2(), # <-- does this need width?
                          linewidth = 0.5 / .pt,
                          width = 0.5) +
            true_line +
            geom_point(position = position_dodge2(width = 0.5),
                       size = 1 / .pt) +
            # geom_hline(yintercept = y_rng[[2]],
            #            colour = "grey", linewidth = 0.5, linetype = "dashed") +
            # scale_colour_manual(breaks = c("uniform", "inverse", "constant"),
            #                     values = c("red", "red", "grey40")) +
            scale_colour_manual(breaks = conv_breaks,
                                values = conv_values) +
            scale_x_discrete(drop = FALSE) +
            scale_y_continuous(limits = ~ range(.x, 0, ymin, ymax)) +
            labs(x = "Scenario",
                 y = "Value",
                 title = html_pars[[par]]) +
            theme_classic() +
            theme(legend.position = "none",
                  plot.title = element_markdown(),
                  axis.text.x = element_text(size = 6,
                                             angle = 45,
                                             hjust = 1))
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
                     rel_heights = c(0.5, pmat$nr))

    if (str_length(alt) > 0) alt <- str_c("-", alt)

    plt_str <- str_glue("{gfx_dir}/{dataset}-all_hpdi{alt}")

    message(str_glue("plotted '{plt_str}'"))

    walk(output, \(op) {
        ggsave(str_glue("{plt_str}.{op}"), plt,
               width = 3 * pmat$nc,
               height = 2 * (pmat$nr + 0.5))
    })

    plt
}

if (FALSE) {
    pars_errorbars("fb-test", 0, "Testing BICI on FB data")
    pars_errorbars("sim-test-inf1", 0, "Validating BICI on Simulated data")
    pars_errorbars("sim-test-inf2", 0, "Validating BICI on Simulated data")
    pars_errorbars("sim-events", 0, "Testing BICI's handling of events")
}
