{
    library(data.table)
    library(purrr)
    library(stringr)
    library(HDInterval)
    library(ggplot2)
    library(ggtext)
    library(ggeasy)
    # library(cowplot)
    
    source("make_parameters.R")
    source("rename_pars.R")
    source("add_h2_to_pars.R")
    source("utils.R")
}


#' Plot posterior distribution
#'
#' Plot a grid of posterior distributions for parameters from a data set.
#'
#' @param dataset name of the data set to be worked on (e.g. `"fb-final"`)
#' @param scen scenario number
#' @param ci credible interval: `"hpdi"` is the highest posterior density
#'   interval. Anything else gives the 95% credible interval
#' @param draw choose between a `"density"` plot or a `"histogram"`
#'
#' @return a list containing various plots and subplots
#' @export
#'
#' @examples
#' `plot_posteriors("fb-final", 1)`
#' `plot_posteriors("fb-final", 1, "hpdi", "histogram")`

plot_posteriors <- function(dataset = "fb-final", scen = 1, rep = 1,
                            ci = "hpdi", # ci95
                            draw = "density" # histogram
) {
    
    # dataset <- "fb-test"; name <- "1-1"; ci <- "hpdi"; draw <- "density"; combine <- TRUE
    # dataset <- "sim-test2"; name <- "1-1"; ci <- "hpdi"; draw <- "density"; combine <- TRUE
    
    name <- str_glue("{scen}-{rep}")
    
    if (rep == 0) name <- scen
    
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
    x[, str_subset(names(x), "state|State|Group") := NULL]
    
    if (x[, all(is.na(map(.SD, sd)))]) {
        message("- All columns are NA!")
        return(NULL)
    }
    
    
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
    
    # Get the list of columns, and find the subset we want to actually plot by
    # discarding unwanted columns like Group effect and shape
    pars <- names(x) |>
        str_subset("^G_|Group|state|State|L_|Prior|Posterior|Number|log", negate = TRUE)
    
    # Remove trial and nested weights if only 1 trial
    if (!str_detect(params$setup, "12")) {
        pars <- if (str_detect(params$setup, "_1")) {
            str_subset(pars, "trial|Tr2|weight2", negate = TRUE)
        } else if (str_detect(params$setup, "_2")) {
            str_subset(pars, "trial|Tr1|weight1", negate = TRUE)
        }
    }
    
    x[, setdiff(names(x), pars) := NULL]
    
    # # Check for values with 0 std. dev. and discard them
    # pars <- x[, names(x)[map(.SD, sd) == 0]]
    # 
    # x[, setdiff(names(x), pars) := NULL]
    
    # Ensure everything is double before melting (in case a col is integer)
    x[, names(x) := map(.SD, as.numeric)]
    
    # Calculate heritability
    pars <- add_h2_to_pars(x, pars)
    h2_priors <- data.table(parameter = c("cov_P_ss", "cov_P_ii", "cov_P_tt",
                                          "h2_ss", "h2_ii", "h2_tt"),
                            type = "Flat", val1 = 0, val2 = 1,
                            true_val = 0.5, use = TRUE)
    priors <- rbind(params$priors, h2_priors)
    
    map(c("ss", "ii", "tt"), \(trait) {
        v2 <- priors[str_detect(parameter, str_c("cov_[GE]_", trait)), val2]
        tv <- priors[str_detect(parameter, str_c("cov_[GE]_", trait)), true_val]
        priors[str_detect(parameter, str_c("cov_P_", trait)),
               `:=`(val2 = sum(v2), true_val = sum(tv))]
    })
    
    
    x2 <- melt(x, measure.vars = pars, variable.name = "parameter")
    
    x3 <- x2[, .(mean = mean(value), median = median(value)), parameter]
    
    # Now convert pars so that it has usable names
    tidy_pars <- rename_pars(pars) |>
        str_replace(" (Tr|Don|Rec)", "\n\\1") |>
        setNames(pars)
    
    plts <- map(pars, \(par) {
        # par <- pars[[1]]
        tidy_param <- tidy_pars[[par]]
        
        xp <- x2[parameter == par, .(value)]
        
        if (sd(xp$value) == 0) {
            # latent_period breaks since sd is 0
            dd <- data.table(x = xp$value[[1]] * c(0.95, 1, 1.05),
                             y = c(0, 1, 0))
        } else {
            dens <- density(xp$value, bw = "SJ", adjust = 1, cut = 0)
            dd <- with(dens, data.table(x, y))
        }
        
        x_mean <- mean(xp$value)
        x_median <- median(xp$value)
        # x_mode <- dd[, x[which.max(y)]] # mode
        
        if (ci == "hpdi") {
            hpdi <- hdi(xp$value, credMass = 0.95)
            lq <- hpdi[["lower"]]
            uq <- hpdi[["upper"]]
        } else {
            lq <- quantile(xp$value, 0.025)
            uq <- quantile(xp$value, 0.975)
        }
        
        par2 <- if (par %in% priors$parameter) {
            par
        } else {
            message(str_glue("Need to fix parameter '{par}'"))
            par |> str_remove_all(c("_Tr[12].*|_Don|_Rec"))
        }
        
        if (priors[parameter == par2, type == "Fixed"]) {
            x_min <- 0
            x_max <- priors[parameter == par2, val1 * 2]
            x_true <- priors[parameter == par2, true_val]
        } else {
            x_min <- priors[parameter == par2, val1]
            x_max <- priors[parameter == par2, val2]
            x_true <- priors[parameter == par2, true_val]
        }
        
        if (draw == "histogram") {
            xbreaks <- seq(x_min, x_max, length.out = 51)
            p <- ggplot() +
                geom_histogram(data = xp, aes(x = value),
                               breaks = xbreaks, fill = "royalblue", colour = NA, alpha = 1) +
                geom_histogram(data = subset(xp, value > lq & value < uq), aes(x = value),
                               breaks = xbreaks, fill = "tomato", colour = NA, alpha = 1)
        } else {
            dd[, zone := factor(fcase(x < lq, "left",
                                      x > uq, "right",
                                      default = "mid"),
                                levels = c("left", "mid", "right"))]
            p <- ggplot(dd) +
                # geom_line(aes(x = x, y = y)) +
                geom_area(aes(x = x, y = y, fill = zone)) +
                scale_fill_manual(values = c("royalblue", "tomato", "royalblue"),
                                  labels = NULL, drop = FALSE)
        }
        
        if (str_starts(dataset, "sim-")) {
            # If we simulated then we know the true underlying parameter
            p <- p + geom_vline(xintercept = x_true,
                                colour = "green")
        }
        
        p <- p + geom_vline(xintercept = x_median,
                            colour = "blue") +
            coord_cartesian(xlim = c(x_min, x_max)) +
            labs(x = "value",
                 y = "density",
                 title = tidy_param) +
            theme_classic() +
            theme(legend.position = "none") +
            easy_title_size(size = 12) +
            easy_remove_y_axis()
        
        p
    }) |> setNames(pars)
    
    title_plt <- ggplot() +
        labs(title = title_plot_str) +
        theme_classic()
    
    # plt_pars <- plot_grid(
    #     title_plt,
    #     plot_grid(plotlist = plts), #ncol = 3),
    #     ncol = 1, rel_heights = c(0.06, 1))
    
    list(title_plt = title_plt,
         plts = plts,
         # pars = plt_pars,
         dataset = dataset,
         scenario = scen)
}
