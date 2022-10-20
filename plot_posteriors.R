{
    library(data.table)
    library(glue)
    library(HDInterval)
    library(ggplot2)
    library(ggtext)
    library(cowplot)
    
    source("rename_pars.R")
}

plot_posteriors <- function(data_set = "fb-mpi",
                            scenario = 1,
                            ci = "hpdi",
                            draw_hist = TRUE) {
    # data_set <- "fb-mpi"; scenario <- 1; ci <- "hpdi"; draw_hist = TRUE
    
    data_dir <- glue("data/{data_set}")
    res_dir <- glue("results/{data_set}")
    
    x <- fread(glue("{data_dir}/scen-{scenario}-1_out/trace_combine.txt"))
    
    res_file <- glue("{res_dir}/scen-{scenario}-1.RData")
    if (file.exists(res_file)) {
        load(res_file)
        title_plot_str <- glue("Scenario {params$label}, {params$description}")
        burn_prop <- with(params, burnin / nsample)
    } else {
        title_plot_str <- glue("{data_set}, scenario {scenario}")
        burn_prop <- 0.2
    }
    
    pars <- parameter_estimates$parameter
    # discard group effect and shape
    pars <- pars[!grepl("^Group effect", pars)]
    pars <- pars[!grepl("_shape$", pars)]

    x <- x[, ..pars]
    
    pars2 <- rename_pars(pars)
    
    x2 <- melt(x, measure.vars = pars, variable.name = "parameter", value.name = "val")
    
    x3 <- x2[, .(mean = mean(val)), parameter]
    
    plts <- list()
    # plts1 <- list()
    for (i in seq_along(pars)) {
        param <- pars[i]
        param2 <- pars2[i]
        
        xp <- x2[parameter == param, .(val)]
        
        dens <- density(xp$val, bw = "SJ", adjust = 1, cut = 0)
        dd <- with(dens, data.table(x, y))
        
        x_mean <- mean(xp$val)
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
        plts[[param]] <- p + geom_vline(xintercept = x_mode, colour = "red") +
            coord_cartesian(xlim = c(x_min, x_max)) +
            labs(x = "value", y = "density", title = param2) +
            theme(legend.position = "none",
                  axis.ticks.y = element_blank(),
                  axis.text.y = element_blank(),
                  axis.title.y = element_blank())
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


data_set <- "fb-mpi"
scenario <- 3

plts <- list()
for (i in 1:3) {
    plts[[i]] <- plot_posteriors(data_set = data_set, scenario = i, ci = "hpdi")
    
    n <- ceiling(length(plts[[i]]$plts) / 3)
    
    pdf(glue("gfx/{data_set}/{data_set}-posteriors-scen{i}.pdf"), width = 10, height = 6)
    print(plts[[i]]$pars)
    dev.off()
}
