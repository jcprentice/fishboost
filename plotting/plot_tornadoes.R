{
    library(data.table)
    library(purrr)
    library(stringr)
    library(lubridate)
    library(HDInterval)
    library(ggplot2)
    library(ggeasy)
    library(cowplot)
    library(extrafont)
    
    source("rename_pars.R")
    source("add_h2_to_pars.R")
    source("utils.R")
}

plot_tornadoes <- function(dataset = "sim-test",
                           scen = 1,
                           sort_by = "median",
                           use_hpdi = TRUE,
                           combine = TRUE) {
    
    # dataset <- "sim-base-inf"; scen <- 1; sort_by <- TRUE; use_hpdi <- TRUE; combine <- TRUE
    
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
    
    trace_files <- list.files(data_dir,
                              pattern = "trace",
                              full.names = TRUE,
                              recursive = TRUE) |>
        str_subset(str_glue("scen-{scen}-")) |>
        str_subset("extended", negate = TRUE) |>
        str_sort(numeric = TRUE)
    
    if (is_empty(trace_files)) {
        message(str_glue("No trace files available for Scenario {scen}!"))
        return(NULL)
    }
    
    
    params <- readRDS(res_files[[1]])$params
    with(params, message(str_glue("Fitted {model_type} to {dataset} dataset",
                         dataset = if (sim_new_data == "no") "Fishboost" else "simulated")))
    
    # Get list of parameters ----
    parameters <- readRDS(res_files[[1]])$parameter_estimates$parameter |>
        # discard group effect and shape
        str_subset("^Group|^G_", negate = TRUE) |>
        rename_bici_pars()
    
    
    # Combine trace files into X ----
    
    # burnin period is 1:(burnin/thin), need to start from (burnin/thin)+1
    burn_prop <- if ("burnprop" %in% names(params)) {
        params$burnprop
    } else {
        with(params, burnin / nsample)
    }
    
    # Some older runs only had a single trace
    if (!any(str_detect(trace_files, "combine"))) combine <- FALSE
    # burnin already discarded in trace_combine
    if (combine) burn_prop <- 0
    trace_files <- str_subset(trace_files, "combine", negate = !combine)
    
    X <- map(trace_files, \(f) {
        # f <- trace_files[[1]]
        tmp <- fread(f)[seq(burn_prop * .N + 1L, .N)]
        tmp[, str_subset(names(tmp), "state|State|^G|^L_|Prior|phi|Number") := NULL]
        
        # Ensure at least 10% of the samples succeeded
        # min_samples <- if (combine) 800 * params$nchains else 800
        # if (nrow(tmp) < min_samples) {
        #     message(str_glue("Too few samples ({nrow(tmp)}) in trace file '{f}'"))
        #     return(NULL)
        # }
        
        tmpx = data.table(parameter = parameters,
                          mean   = tmp[, map(.SD, mean)]   |> unlist() |> unname(),
                          median = tmp[, map(.SD, median)] |> unlist() |> unname())
        
        CIs <- if (use_hpdi) {
            # 95% Highest Density Posterior Interval
            tmp[, map(.SD, hdi, credMass = 0.95)]
        } else {
            # 95% Credible Interval (with 2.5% tails)
            tmp[, map(.SD, quantile, c(0.025, 0.975))]
        }
        tmpx[, `:=`(min = CIs[1] |> unlist() |> unname(),
                    max = CIs[2] |> unlist() |> unname())]
        tmpx
    }) |>
        rbindlist(idcol = "chain")
    
    
    # Tidy and sort X ----
    if (sort_by %in% c("mean", "median")) {
        setorderv(X, cols = c("parameter", sort_by))
    } else {
        setorder(X, "parameter")
    }
    X[, chain := seq(.N), by = parameter] |>
        setnames("chain", "id", skip_absent = TRUE)
    
    
    # Get true and posterior means ----
    
    Xtab <- data.table(parameter = parameters, xmin = NA_real_, xmax = NA_real_,
                       true_val = NA_real_, est_val = NA_real_,
                       hdi1 = NA_real_, hdi2 = NA_real_, id = last(X$id) + 1L)
    
    iwalk(parameters, \(par, i) {
        # par2 <- str_remove_all(par, "_Tr.*")
        X_mu <- X[parameter == par, mean(mean)]
        X_hdi <-  X[parameter == par, hdi(mean)][c("lower", "upper")]
        set(Xtab, i, c("xmin", "xmax", "true_val", "est_val", "hdi1", "hdi2"),
            c(params$priors[parameter == par, .(val1, val2, true_val)], X_mu, X_hdi))
    })
    
    # Set X-axis limits ----
    Xtab[str_starts(parameter, "cov_"), `:=`(xmin = 0, xmax = pmax(1, xmax))]
    Xtab[str_starts(parameter, "r_"), `:=`(xmin = -1, xmax = 1)]
    Xtab[str_detect(parameter, "beta|period|shape|period"), xmin := 0]
    
    # print(Xtab)
    
    # Prevent drawing green "true" line if using experimental data
    if (str_detect(dataset, "fb")) Xtab[, true_val := NA_real_]
    
    # Make these nice to read when plotted
    params2 <- setNames(rename_pars(parameters), parameters)
    
    
    # Plot parameters ----
    plots <- map(parameters, \(par) {
        par2 <- params2[[par]]
        x1 <- Xtab[parameter == par] |> as.list()
        
        ggplot(data = X[parameter == par]) +
            # credible interval
            geom_segment(aes(x = min, xend = max, y = id, yend = id, colour = par)) +
            # mean
            geom_point(aes(x = mean, y = id),
                       colour = "red", size = 1) +
            # true parameter value
            geom_vline(xintercept = x1$true_val, colour = "green") +
            geom_vline(xintercept = x1$est_val, colour = "blue") +
            geom_vline(xintercept = x1$hdi1, colour = "blue",
                       linetype = "dashed", linewidth = 0.5) +
            geom_vline(xintercept = x1$hdi2, colour = "blue",
                       linetype = "dashed", linewidth = 0.5) +
            coord_cartesian(xlim = c(x1$xmin, x1$xmax)) +
            expand_limits(x = c(0, x1$true_val)) +
            labs(title = par2, x = "Value", y = "id") +
            theme_classic() +
            theme(plot.title = element_text(size = 10)) +
            easy_remove_legend() +
            easy_remove_axes("y")
    }) |>
        setNames(parameters)
    

    label <- params$label %||% scen
    description <- params$description |>
        str_remove(", (coverage|convergence)")

    title_plt <- ggplot() +
        labs(title = str_glue("{dataset} / {label}: {description}")) +
        theme_classic()
    
    
    # Combine plots into grid ----
    
    plt_pars <- plot_grid(
        title_plt,
        plot_grid(plotlist = plots), # ncol = 3?
        ncol = 1, rel_heights = c(0.05, 1))
    
    # return
    list(X = X,
         title_plt = title_plt,
         plots = plots,
         pars = plt_pars)
}


# https://ben-williams.github.io/updated_ggplot_figures.html

# font_import() only do this one time - it takes a while
# loadfonts(device = "pdf")


# For the WCGALP abstract:
#
# theme_set(theme_bw(base_size = 12, base_family = "Times New Roman") +
#               theme(panel.grid.major = element_blank(),
#                     panel.grid.minor = element_blank()))

