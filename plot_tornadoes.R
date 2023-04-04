{
    library(data.table)
    library(glue)
    library(gtools)
    library(lubridate)
    library(HDInterval)
    library(ggplot2)
    library(cowplot)
    library(extrafont)
    
    source("rename_pars.R")
    source("utils.R")
}

plot_tornadoes <- function(data_set = "sim-linked",
                           scenario = 1,
                           sort_CIs = TRUE,
                           use_hpdi = TRUE,
                           combine = TRUE,
                           show_plots = FALSE) {
    
    # data_set = "sim-censored-1"; scenario = 1; sort_CIs = TRUE; use_hpdi = TRUE; combine = FALSE; show_plots = FALSE
    
    
    # slightly different behaviour if using actual data since we don't know the true means
    fb <- grepl("fb", data_set)
    mpi <- grepl("mpi", data_set)
    
    # Directory to save all files to
    gfx_dir <- glue("gfx/{data_set}")
    if (!dir.exists(gfx_dir)) {
        dir.create(gfx_dir)
    }
    
    message(glue("Generating plots for {data_set}, scenario {scenario} ..."))
    
    # get list of files to scan
    data_dir  <- glue("data/{data_set}")
    res_dir   <- glue("results/{data_set}")
    res_files <- list.files(res_dir, pattern = glue("scen-{scenario}-"))
    
    trace_files <- list.files(data_dir, pattern = "trace", recursive = TRUE, include.dirs = TRUE)
    trace_files <- trace_files[grepl(glue("scen-{scenario}-"), trace_files)]
    trace_files <- mixedsort(trace_files, decreasing = TRUE)
    
    if (length(res_files) == 0L) {
        message(glue("No files available for Scenario {scenario}!"))
        return(NULL)
    } else {
        message(glue("Found {length(res_files)} file(s)"))
    }
    
    run_times <- numeric(length(res_files))
    
    for (i in seq_along(res_files)) {
        # Loads: params, pop, parameter_estimates, estimated_BVs, pred_accs,
        # ranks, time_taken
        f <- glue("{res_dir}/{res_files[i]}")
        load(f)
        
        if (i == 1) {
            message(glue("Fitted {mt} to {ds} dataset",
                         mt = params$model_type,
                         ds = if (params$use_fb_data) "Fishboost" else "simulated"))
        }
        
        run_times[i] <- time_taken[["elapsed"]]
    }
    
    # Check how much time we needed here
    message(glue("Runs took about {a} (longest {b})",
                 a = seconds_to_period(ceiling(mean(run_times))),
                 b = seconds_to_period(ceiling(max(run_times)))))
    
    
    # Get list of parameters ----
    parameters <- parameter_estimates$parameter
    # discard group effect and shape
    parameters <- parameters[!grepl("^Group effect", parameters)]
    parameters <- parameters[!grepl("_shape$", parameters)]
    
    
    # Combine trace files into X ----
    X <- data.table()
    
    # burnin period is 1:(burnin/thin), need to start from (burnin/thin)+1
    burn_prop <- if (exists("params")) with(params, burnin / nsample) else 0.2
    burn_prop <- clamp(burn_prop, 0, 1)
    
    if (params$sire_version == "sire22" && combine) {
        trace_files <- trace_files[grepl("combine", trace_files)]
        # burnin was already discarded
        burn_prop <- 0
    } else {
        trace_files <- trace_files[!grepl("combine", trace_files)]
    }
    
    for (tf in trace_files) {
        tmp <- fread(glue("{data_dir}/{tf}"))[(ceiling(burn_prop * .N) + 1):.N]
        
        # Ensure at least 10% of the samples succeeded
        if (nrow(tmp) < with(params, 0.1 * (nsample - burnin) / thin)) {
            message(glue("Too few samples ({nrow(tmp)}) in trace file '{tf}'"))
            next
        }
        
        tmpx = data.table(parameter = parameters)
        
        for (param in parameters) {
            if (use_hpdi) {
                # 95% Highest Density Posterior Interval
                hpdi = tmp[, hdi(get(param), credMass = 0.95)]
                tmpx[parameter == param, `:=`(min  = hpdi[["lower"]],
                                              mean = tmp[, mean(get(param))],
                                              max  = hpdi[["upper"]])]
            } else {
                # 95% Credible Interval (with 2.5% tails)
                tmpx[parameter == param, `:=`(min  = tmp[, quantile(get(param), .025)],
                                              mean = tmp[, mean(get(param))],
                                              max  = tmp[, quantile(get(param), .975)])]
            }
        }
        
        X <- rbind(X, tmpx)
    }
    
    
    
    # Tidy and sort X ----
    if (sort_CIs) {
        setkey(X, "parameter", "mean")
    } else {
        setkey(X, "parameter")
    }
    X[, id := seq(.N), by = parameter]
    
    
    # Get true and posterior means ----
    
    par_means <- params$priors[parameter %in% parameters, .(parameter, true_val)]
    par_means <- merge(par_means,
                       X[parameter %in% parameters, .(est_val = mean(mean)), by = parameter],
                       by = "parameter")
    
    # Set X-axis limits ----
    # Xrng <- params$priors[parameter %in% parameters, .(parameter, xmin = val1, xmax = val2)]
    Xrng <- X[, .(xmin = min(min), xmax = max(max)), by = parameter]
    
    # Xrng[startsWith(parameter, "r_"),
    #      `:=`(xmin = min(-1, xmin), xmax = max(1, xmax)),
    #      by = parameter]
    
    Xrng[startsWith(parameter, "cov_"),
         `:=`(xmin = 0, xmax = max(1, xmax)),
         by = parameter]
    
    Xrng[startsWith(parameter, "r_"),
         `:=`(xmin = -1, xmax = +1),
         by = parameter]
    
    
    Xrng[parameter %in% c("beta", "recovery_period", "gamma_shape",
                          "latency_period", "eta_shape",
                          "detection_period", "rho_shape"),
         xmin := 0]
    # xmax = max(xmax, par_means[parameter, true_val], something like that?
    Xrng[parameter == "σ_G", `:=`(xmin = 0, xmax = xmax)]
    
    # Xrng[parameter == "beta", xmax:= 0.5]
    
    
    # Prevent drawing green "true" line if using real data
    if (fb) par_means[, true_val := NA_real_]
    
    print(Xrng)
    
    # Make these nice to read when plotted
    params2 <- rename_pars(parameters)
    
    
    # Plot parameters ----
    plots <- list()
    
    for (i in seq_along(parameters)) {
        param <- parameters[i]
        param2 <- params2[i]
        df <- X[parameter == param]
        
        pm <- par_means[parameter == param]
        pm[, id := nrow(df) + 1]
        pm <- as.list(pm[1])
        
        plots[[param]] <- ggplot(data = df) +
            # credible interval
            geom_segment(aes(x = min, xend = max, y = id, yend = id, colour = param)) +
            # mean
            geom_point(aes(x = mean, y = id, colour = param)) +
            # true parameter value
            geom_vline(xintercept = pm$true_val, colour = "green") +
            geom_vline(xintercept = pm$est_val, colour = "blue") +
            coord_cartesian(xlim = c(min(Xrng[param, xmin], 0),
                                     max(Xrng[param, xmax], 0))) +
            labs(title = param2, x = "Value", y = "id") +
            theme(legend.position = "none",
                  axis.ticks.y = element_blank(),
                  axis.text.y = element_blank(),
                  axis.title.y = element_blank())
    }
    
    if (is.null(params$label)) {
        params$label <- scenario
    }
    title_plt <- ggplot() + labs(title = glue("Scenario {params$label}, {params$description}"))
    
    
    # Combine plots into grid ----
    
    plt_pars <- plot_grid(
        title_plt,
        plot_grid(plotlist = plots), # ncol = 3?
        ncol = 1, rel_heights = c(0.06, 1))
    
    pars_str <- glue("{gfx_dir}/pars-scen{scenario}.png")
    message("Plotting Parameters to ", pars_str)
    ggsave(pars_str, plt_pars, width = 2000, height = 1000, units = "px")
    
    # Save data ----
    save(X, plt_pars,
         file = glue("meta/violins-{data_set}-{scenario}.RData"))
    
    if (show_plots) {
        print(plt_rank)
        print(plt_pars)
    }
    
    # return
    list(X = X,
         title_plt = title_plt,
         plots = plots,
         pars = plt_pars)
}


# run tests ----

# plts1 <- plot_tornadoes(data_set = "sim", scenario = 1, show_plots = FALSE)
# plts1 <- plot_tornadoes(data_set = "fb_1_12_Dx_linked", scenario = 17, show_plots = FALSE)

# data_set <- "fb-mpi"; combine <- FALSE
data_set <- "sim-censored-1"; combine <- FALSE
# data_set <- "sim-Gsi_cov_Da-1-mpi"; combine = FALSE
# data_set = "sim-donor_links1-2-mpi"; combine = TRUE

if (TRUE) {
    plts_sim <- list()
    for (i in 1:3) {
        # combine: 1=FALSE, 2=TRUE
        plts_sim[[i]] <- plot_tornadoes(data_set = data_set, combine = combine, scenario = i)
        
        # p_ht <- ceiling(sqrt(length(plts_sim[[i]]$plots))) * 5/3
        
        ggsave(glue("gfx/{data_set}/{data_set}-pars-scen{i}.pdf"),
               plts_sim[[i]]$pars,
               width = 10, height = 6)
    }
}





# https://ben-williams.github.io/updated_ggplot_figures.html

# font_import() only do this one time - it takes a while
# loadfonts(device = "pdf")


# For the WCGALP abstract:
#
# theme_set(theme_bw(base_size = 12, base_family = "Times New Roman") +
#               theme(panel.grid.major = element_blank(),
#                     panel.grid.minor = element_blank()))
