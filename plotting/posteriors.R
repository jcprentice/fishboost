{
    library(data.table)
    library(glue)
    library(HDInterval)
    library(ggplot2)
    library(ggtext)
    library(cowplot)
    
    source("rename_pars.R")
    source("add_h2.R")
}

plot_posteriors <- function(
        data_set = "fb-parasites",
        scenario = 1,
        ci = "hpdi",
        draw_hist = TRUE
) {
    
    # data_set <- "fb-parasites4"; scenario <- 1; ci <- "hpdi"; draw_hist = TRUE
    
    # Slight nudge to ensure that LP is the same type as everything else
    x <- fread(glue("data/{data_set}/scen-{scenario}-1_out/trace_combine.txt"))
    
    res_file <- glue("results/{data_set}/scen-{scenario}-1.RData")
    if (file.exists(res_file)) {
        load(res_file)
        title_plot_str <- glue("Scenario {params$label}, {params$description}")
        burn_prop <- with(params, burnin / nsample)
    } else {
        title_plot_str <- glue("{data_set}, scenario {scenario}")
        burn_prop <- 0.2
    }
    
    # Get the list of columns, and find the subset we want to actually plot
    pars <- names(x)
    
    # Discard unwanted columns like group effect and shape
    pars <- pars[!grepl("^Group effect", pars)]
    pars <- pars[!grepl("_shape$", pars)]
    pars <- pars[!grepl("state", pars)]
    
    x <- x[, names(x)[!names(x) %in% pars] := NULL]
    
    # Check for values with 0 std. dev. and discard them
    y <- x[, lapply(.SD, sd)][, .(.SD != 0), .SDcols = pars]
    pars <- names(y)[y == TRUE]
    
    x[, names(x)[!names(x) %in% pars] := NULL]
    
    # Calculate heritability
    pars <- add_h2(x, pars)
    h2_priors <- data.table(parameter = c("h2_ss", "h2_ii", "h2_rr"),
                            type = "Flat", val1 = 0, val2 = 1, true_val = 0.5, use = TRUE)
    params$priors <- rbind(params$priors, h2_priors)
    
    
    x2 <- melt(x, measure.vars = pars, variable.name = "parameter", value.name = "val")
    
    x3 <- x2[, .(mean = mean(val)), parameter]
    
    # Now convert pars so that it has usable names
    tidy_pars <- rename_pars(pars)
    names(tidy_pars) <- pars
    
    plts <- list()
    # plts1 <- list()
    for (param in pars) {
        tidy_param <- tidy_pars[[param]]
        
        xp <- x2[parameter == param, .(val)]
        
        dens <- density(xp$val, bw = "SJ", adjust = 1, cut = 0)
        dd <- with(dens, data.table(x, y))
        
        x_mean <- mean(xp$val)
        x_median <- median(xp$val)
        x_mode <- dd[, x[which.max(y)]] # mode
        
        if (ci == "hpdi") {
            hpdi <- hdi(xp$val, credMass = 0.95)
            lq <- hpdi[["lower"]]
            uq <- hpdi[["upper"]]
        } else {
            lq <- quantile(xp$val, 0.025)
            uq <- quantile(xp$val, 0.975)
        }
        
        x_min <- params$priors[parameter == param, val1]
        x_max <- params$priors[parameter == param, val2]
        x_true <- params$priors[parameter == param, true_val]
        
        
        if (draw_hist) {
            xbreaks <- seq(x_min, x_max, length.out = 51)
            p <- ggplot() +
                geom_histogram(data = xp, aes(x = val),
                               breaks = xbreaks, fill = "royalblue", colour = NA, alpha = 1) +
                geom_histogram(data = subset(xp, val > lq & val < uq), aes(x = val),
                               breaks = xbreaks, fill = "tomato", colour = NA, alpha = 1)
        } else {
            p <- ggplot() +
                geom_line(data = dd, aes(x = x, y = y)) +
                geom_ribbon(data = subset(dd, lq < x & x < uq),
                            aes(x = x, ymin = 0, ymax = y),
                            fill = "red", colour = NA, alpha = 0.5)
        }
        if (startsWith(data_set, "fb")) {
            p <- p +
                geom_vline(xintercept = x_mean, colour = "red") +
                geom_vline(xintercept = x_median, color = "green") +
                geom_vline(xintercept = x_mode, colour = "blue")
        } else {
            p <- p +
                geom_vline(xintercept = x_mean, colour = "blue") +
                geom_vline(xintercept = x_true, colour = "green")
        }
        p <- p +
            coord_cartesian(xlim = c(x_min, x_max)) +
            labs(x = "value", y = "density", title = tidy_param) +
            theme(legend.position = "none",
                  axis.ticks.y = element_blank(),
                  axis.text.y = element_blank(),
                  axis.title.y = element_blank())
        
        plts[[param]] <- p
    }
    
    title_plt <- ggplot() + labs(title = glue("Scenario {params$label}, {params$description}"))
    
    plt_pars <- plot_grid(
        title_plt,
        plot_grid(plotlist = plts), #ncol = 3),
        ncol = 1, rel_heights = c(0.06, 1))
    
    list(title_plt = title_plt,
         plts = plts,
         pars = plt_pars)
}

