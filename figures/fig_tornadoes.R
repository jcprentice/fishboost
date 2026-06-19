{
    library(data.table)
    library(purrr)
    library(stringr)
    library(HDInterval)
    library(ggplot2)
    library(ggtext)
    # library(ggeasy)
    library(cowplot)

    source("rename_pars.R")
    source("figures/theme_jamie.R")
}

fig_tornadoes <- function(dataset = "sim-test-inf1", scen = 1) {
    if (FALSE) {
        dataset <- "sim-test-inf1"
        scen <- 2
    }

    message(str_glue("Generating tornado plots for {dataset} / s{scen} ..."))

    # Directories and list of files to scan
    data_dir    <- str_glue("datasets/{dataset}/data")
    results_dir <- str_glue("datasets/{dataset}/results")
    res_files   <- list.files(results_dir,
                              pattern = str_glue("scen-{scen}-"),
                              full.names = TRUE) |>
        str_sort(numeric = TRUE)

    nres <- length(res_files)
    if (nres == 0) {
        message(str_glue("No results files available for Scenario {scen}!"))
        return(NULL)
    }
    message(str_glue("Found {nres} results file{s}",
                     s = if (nres != 1) "s" else ""))

    params <- readRDS(res_files[[1]])$params

    # Get list of parameters ----

    pars <- readRDS(res_files[[1]])$parameter_estimates$parameter |>
        str_subset("^G", negate = TRUE)

    # May need to rename weights if only 1 trial
    trials <- str_extract(pars, "beta_Tr(.)", group = 1) |>
        discard(is.na) |> unique() |> as.integer() |> sort()

    # Combine trace files into X ----

    trace_files <- list.files(data_dir,
                              pattern = "trace_combine",
                              full.names = TRUE,
                              recursive = TRUE) |>
        str_subset(str_glue("scen-{scen}-")) |>
        str_subset("extended", negate = TRUE) |>
        str_sort(numeric = TRUE)

    if (is_empty(trace_files)) {
        message(str_glue("No trace files available for Scenario {scen}!"))
        return(NULL)
    }

    X <- map(trace_files, \(f) {
        if (FALSE) {
            f <- trace_files[[1]]
        }
        tmp <- fread(f)[, .SD, .SDcols = !patterns("state|^G")]

        tmpx = data.table(parameter = names(tmp),
                          mean   = tmp[, map(.SD, mean)]   |> unlist(),
                          median = tmp[, map(.SD, median)] |> unlist())

        # 95% Highest Density Posterior Interval
        CIs <- tmp[, map(.SD, hdi, credMass = 0.95)]
        tmpx[, `:=`(ci_min = CIs[1] |> unlist(),
                    ci_max = CIs[2] |> unlist())]
        tmpx
    }) |> rbindlist()


    # Tidy and sort X ----

    setorder(X, parameter, median)
    X[, id := seq_len(.N), parameter]

    # Get true and posterior means ----

    Xtab <- data.table(parameter = pars, xmin = NA_real_, xmax = NA_real_,
                       true_val = NA_real_, est_val = NA_real_,
                       hdi1 = NA_real_, hdi2 = NA_real_)

    iwalk(pars, \(par, i) {
        # html_par <- str_remove_all(par, "_Tr.*")
        X_mu <- X[parameter == par, mean(mean)]
        X_hdi <-  X[parameter == par, hdi(mean)][c("lower", "upper")]
        set(Xtab, i, c("xmin", "xmax", "true_val", "est_val", "hdi1", "hdi2"),
            c(params$priors[parameter == par, .(val1, val2, true_val)], X_mu, X_hdi))
    })

    rename_weights <- function(x, n) {
        str_replace(x, "weight_", str_c("weight", n, "_"))
    }

    if (!identical(trials, c(1L, 2L))) {
        pars <- rename_weights(pars, trials)
        X[, parameter := rename_weights(parameter, trials)]
        Xtab[, parameter := rename_weights(parameter, trials)]
    }

    # Set X-axis limits ----
    Xtab[str_starts(parameter, "cov_"), `:=`(xmin = 0, xmax = pmax(1, xmax))]
    Xtab[str_starts(parameter, "r_"), `:=`(xmin = -1, xmax = 1)]
    Xtab[str_detect(parameter, "beta|period|shape|period"), xmin := 0]

    # print(Xtab)

    # Prevent drawing green "true" line if using experimental data
    if (str_detect(dataset, "fb")) Xtab[, true_val := NA_real_]

    # Make these nice to read when plotted
    html_pars <- setNames(html_names(pars), pars)

    l2p <- 1 / ggplot2::.pt

    # Plot parameters ----

    plts <- map(pars, \(par) {
        if (FALSE) {
            par <- pars[[1]]
        }
        html_par <- html_pars[[par]]
        x1 <- Xtab[parameter == par, .SD, .SDcols = -1] |> unlist()

        ggplot(data = X[parameter == par]) +
            geom_segment(aes(x = ci_min, xend = ci_max, y = id, yend = id,
                             colour = par)) +
            geom_point(aes(x = mean, y = id),
                       colour = "red", size = 1 * l2p) +
            geom_vline(xintercept = x1[c("true_val", "est_val", "hdi1", "hdi2")],
                       colour = c("green", "blue", "blue", "blue"),
                       linewidth = c(1, 1, 0.5, 0.5) * l2p,
                       linetype = c("solid", "solid", "dashed", "dashed")) +
            scale_x_continuous(limits = ~ range(.x, 0, x1[c("true_val", "xmin", "xmax")])) +
            labs(x = "Value",
                 y = NULL,
                 title = html_par) +
            theme_jamie() +
            ggeasy::easy_remove_axes("y")
    }) |>
        setNames(pars)

    plt_names <- c(
        "cov_G_ss",   "cov_G_ii",   "cov_G_tt",   "r_G_si",    "r_G_st",    "r_G_it",
        "cov_E_ss",   "cov_E_ii",   "cov_E_tt",   "r_E_si",    "r_E_st",    "r_E_it",
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

    plt_names <- intersect(plt_names, pars)

    plt <- plot_grid(plotlist = plts[plt_names],
                     ncol = 6, align = "hv")
    plt

    plt_str <- str_glue("gfx/{dataset}-s{scen}-tornados")
    ggsave(str_glue(plt_str, ".pdf"), plt,
           width = 18.3, height = 17, units = "cm")
    ggsave(str_glue(plt_str, ".png"), plt,
           width = 18.3, height = 17, units = "cm", dpi = "print")

    plt
}

if (FALSE) {
    plt <- fig_tornadoes("sim-test-inf1", 2)
    plt <- fig_tornadoes("sim-test-inf1", 4)
    plt
}

