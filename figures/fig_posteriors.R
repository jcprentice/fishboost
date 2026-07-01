{
    library(data.table)
    library(purrr)
    library(stringr)
    library(HDInterval)
    library(ggplot2)
    library(cowplot)
    library(ggtext)
    library(ggeasy)
    # library(cowplot)

    source("make_parameters.R")
    source("rename_pars.R")
    source("add_h2_to_pars.R")
    source("utils.R")
    source("figures/theme_natcom.R")
}


#' Plot posterior distribution
#'
#' Plot a grid of posterior distributions for parameters from a data set.
#'
#' @param dataset name of the data set to be worked on (e.g. `"fb-final"`)
#' @param scen scenario number
#' @param draw choose between a `"density"` plot or a `"histogram"`
#'
#' @return a list containing various plots and subplots
#' @export
#'
#' @examples
#' `plot_posteriors("fb-final", 1)`
#' `plot_posteriors("fb-final", 1, "hpdi", "histogram")`

fig_posteriors <- function(dataset = "fb-final", scen = 1, rep = 1) {

    if (FALSE) {
        dataset <- "fb-test"; scen <- 7; rep <- 1
        dataset <- "sim-test2"; scen <- 1; rep <- 1
    }

    name <- if (rep == 0) scen else c(str_glue("{scen}-{rep}"))

    message(str_glue("Generating posterior plots for '{dataset} / s{name}' ..."))

    data_dir    <- str_glue("datasets/{dataset}/data")
    results_dir <- str_glue("datasets/{dataset}/results")

    # Slight nudge to ensure that LP is the same type as everything else
    trace_files <- list.files(data_dir, "trace_combine.tsv",
                        recursive = TRUE, full.names = TRUE) |>
        str_subset(str_glue("scen-{scen}-")) |>
        str_sort(numeric = TRUE)

    if (rep != 0) {
        trace_files <- str_subset(trace_files, str_glue("-{rep}-out"))
    }

    if (is_empty(trace_files)) {
        message("- No trace files found")
        return(NULL)
    }

    x <- map(trace_files, fread) |>
        rbindlist()
    x[, str_subset(names(.SD), "state|State|Group|^G_") := NULL]

    if (x[, all(is.na(map(.SD, sd)))]) {
        message("- All columns are NA!")
        return(NULL)
    }

    trials <- str_extract(names(x), "beta_Tr(.)", group = 1) |>
        discard(is.na) |> unique() |> as.integer() |> sort()


    if (any(str_detect(names(x), "weight_"))) {
        setnames(x, \(s) {
            str_replace_all(s, "weight_", str_glue("weight{trials[[1]]}_"))
        })
    }

    pars <- names(x)

    res_file <- list.files(results_dir, full.names = TRUE) |>
        str_subset(str_glue("scen-{scen}-")) |>
        str_sort(numeric = TRUE) |>
        pluck(1)

    if (!is.null(res_file)) {
        params <- readRDS(res_file)$params
        title_plot_str <- str_glue("{dataset}/s{name}, {params$description}") |>
            str_remove(", (convergence|coverage)")
        if ("sim_new_data" %notin% names(params)) {
            params$sim_new_data <- if (str_detect(dataset, "sim")) "r" else "no"
        }
    } else {
        message("- No results file found, attempting to proceed without one.")
        rf <- str_glue("param_sets/{dataset}.rds")
        if (file.exists(rf)) {
            res <- readRDS(rf)
            protocol <- res$protocol[scenario == scen & replicate == max(rep, 1)] |>
                as.list()

            params <- protocol[c("setup", "description")]
            title_plot_str <- str_glue("{dataset}/{name}, {params$description}") |>
                str_remove(", (convergence|coverage)")
        } else {
            params <- make_parameters(sim_new_data = if(str_detect(dataset, "sim")) "r" else "no")
            title_plot_str <- str_glue("{dataset}/s{name}")
        }
    }

    # # Check for values with 0 std. dev. and discard them
    # pars <- x[, names(x)[map(.SD, sd) == 0]]
    #
    # x[, setdiff(names(x), pars) := NULL]

    # Ensure everything is double before melting (in case a col is integer)
    x[, names(.SD) := map(.SD, as.numeric)]

    # Calculate heritability
    if (!any(str_detect(pars, "h2")) && any(str_detect(pars, "_G_"))) {
        message(str_glue("- No heritability found, suggest running\n",
                         "> rebuild_bici_posteriors(\"{x}\", \"{y}\")",
                         x = dataset, y = str_glue("scen-{name}")))
        pars <- add_h2_to_pars(x, pars)
    }

    priors <- rbind(params$priors,
                    data.table(parameter = c("cov_P_ss", "cov_P_ii", "cov_P_tt",
                                             "h2_ss", "h2_ii", "h2_tt"),
                               type = "default", val1 = 0, val2 = 1,
                               true_val = 0.5, use = TRUE))

    walk(c("ss", "ii", "tt"), \(trait) {
        v2 <- priors[str_detect(parameter, str_c("cov_[GE]_", trait)), val2]
        tv <- priors[str_detect(parameter, str_c("cov_[GE]_", trait)), true_val]
        priors[str_detect(parameter, str_c("cov_P_", trait)),
               `:=`(val2 = sum(v2), true_val = sum(tv))]
    })

    priors[str_starts(parameter, "r_"), `:=`(val1 = -1, val2 = 1)]
    priors[type == "inverse",
           val1 := fcase(str_detect(parameter, "beta"), pmax(0.1, val1),
                         str_detect(parameter, "[LDR]P"), pmax(1, val1),
                         default = pmax(0.01, val1))]


    x2 <- melt(x, measure.vars = pars, variable.name = "parameter")

    x3 <- x2[, .(mean = mean(value), median = median(value)), parameter]

    # Now convert pars so that it has usable names
    html_pars <- html_names(pars) |>
        # str_replace(" (\\(.+\\))", "\n\\1") |>
        setNames(pars)

    plts <- map(pars, \(par) {
        if (FALSE) {
            par <- pars[[1]]
        }
        html_par <- html_pars[[par]]

        xp <- x2[parameter == par, .(value)]

        if (sd(xp$value) == 0) {
            # LP breaks since sd is 0
            dd <- data.table(x = xp$value[[1]] * c(0.95, 1, 1.05),
                             y = c(0, 1, 0))
        } else {
            dens <- density(xp$value, bw = "SJ", adjust = 1, cut = 0)
            dd <- with(dens, data.table(x, y))
        }

        x_mean <- mean(xp$value)
        x_median <- median(xp$value)
        # x_mode <- dd[, x[which.max(y)]] # mode

        hpdi <- hdi(xp$value, credMass = 0.95)
        lq <- hpdi[["lower"]]
        uq <- hpdi[["upper"]]

        par2 <- if (par %in% priors$parameter) {
            par
        } else {
            message(str_glue("Need to fix parameter '{par}'"))
            par |> str_remove_all(c("_Tr[12].*|_Don|_Rec"))
        }

        prior_type <- priors[parameter == par2, type]

        if (prior_type == "Fixed") {
            x_min  <- 0
            x_max  <- priors[parameter == par2, val1 * 2]
            x_true <- priors[parameter == par2, true_val]
        } else {
            x_min  <- priors[parameter == par2, val1]
            x_max  <- priors[parameter == par2, val2]
            x_true <- priors[parameter == par2, true_val]
        }

        dd[, zone := factor(fcase(x < lq, "left",
                                  x > uq, "right",
                                  default = "mid"),
                            levels = c("left", "mid", "right"))]

        prior <- data.table(x = seq(x_min, x_max, length.out = 101))
        if (prior_type == "uniform") {
           prior[, y := 1 / (x_max - x_min)]
        } else if (prior_type == "inverse") {
           prior[, y := pmin(1 / x)]
        } else if (prior_type == "default") {
           if (str_detect(par, "cov_")) {
               prior[, y := dnorm(x, 0, 2) * 2]
           } else if (str_detect(par, "r_")) {
               prior[, y := dbeta((x + 1) / 2, 1.2, 1.2) / 2]
           } else if (str_detect(par, "h2_")) {
               prior[, y := NA]
           }
        } else {
            message(str_glue("Define prior for parameter {par} with ",
                             "prior_type = '{prior_type}' "))
        }
        prior[, y := y / sum(c(0, diff(x)) * y, na.rm = TRUE)]

        l2p <- 1 / ggplot2::.pt

        true_val <- if (str_starts(dataset, "sim-")) {
            geom_vline(xintercept = x_true,
                       colour = "green",
                       linewidth = 1 * l2p)
        }
        post_median <-  geom_vline(xintercept = x_median,
                                   colour = "blue",
                                   linewidth = 1 * l2p)

        p <- ggplot(dd) +
            # geom_line(aes(x = x, y = y)) +
            geom_area(aes(x = x, y = y, fill = zone)) +
            true_val +
            post_median +
            geom_line(data = prior, aes(x = x, y = y),
                      colour = "black",
                      linetype = "dashed") +
            scale_fill_manual(values = c("royalblue", "tomato", "royalblue"),
                              labels = NULL, drop = FALSE) +
            coord_cartesian(xlim = c(x_min, x_max)) +
            labs(x = "value",
                 y = "density",
                 title = html_par) +
            theme_natcom() +
            theme(legend.position = "none") +
            easy_remove_y_axis()
        p
    }) |> setNames(pars)

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

    plot_str <- str_glue("gfx/{dataset}-s{scen}-posteriors")
    ggsave(str_glue("{plot_str}.pdf"), plt,
           width = 18.3, height = ht, units = "cm")
    ggsave(str_glue("{plot_str}.png"), plt,
           width = 18.3, height = ht, units = "cm", dpi = "print")

    plt
}

fig_posteriors("fb-test", 7, 1)
