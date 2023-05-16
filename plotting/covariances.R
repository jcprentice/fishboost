library(data.table)
library(glue)
library(HDInterval)
library(ggplot2)
library(ggtext)
library(cowplot)


plot_covariances <- function(data_set = "fb-parasites",
                             scenario = 1,
                             ci = "hpdi",
                             replicate = 0) {
    
    # data_set <- "fb-parasites"; scenario <- 1; ci <- "hpdi"; replicate <- 0
    
    res_files <- dir(glue("results/{data_set}"),
                     pattern = glue("scen-{scenario}-1.RData"),
                     full.names = TRUE)
    if (length(res_files) > 0) {
        load(res_files[1])
        title_plot_str <- glue("Scenario {params$label}, {params$description}")
        burn_prop <- with(params, burnin / nsample)
    } else {
        title_plot_str <- glue("{data_set}, scenario {scenario}")
        burn_prop <- 0.2
    }
    
    # Names as in trace
    pars <- c("cov_G_ss", "cov_G_ii", "cov_G_rr", "r_G_si", "r_G_sr", "r_G_ir")
    pars <- c(pars, sub("G", "E", pars), "sigma")
    
    # Names as wanted on figure
    pars2 <- c("var<sub>G</sub>(g)", "var<sub>G</sub>(f)", "var<sub>G</sub>( r)",
               "cor<sub>G</sub>(g,f)", "cor<sub>G</sub>(g,r)", "cor<sub>G</sub>(f,r)")
    pars2 <- c(pars2, sub("G", "E", pars2), "sigma")
    
    
    trace_file <- glue("data/{data_set}/scen-{scenario}-1_out/trace_combine.txt")
    x <- fread(trace_file)[, ..pars]
    
    
    x[, `:=`(
        h2_G_ss = cov_G_ss / (cov_G_ss + cov_E_ss),
        h2_G_ii = cov_G_ii / (cov_G_ii + cov_E_ii),
        h2_G_rr = cov_G_rr / (cov_G_rr + cov_E_rr)
    )]
    
    x2 <- melt(x, measure.vars = pars, variable.name = "parameter", value.name = "val")
    
    x3 <- x2[, .(median = median(val)), parameter]
    
    plts <- list()
    # plts1 <- list()
    for (i in seq_along(pars)) {
        param <- pars[i]
        param2 <- pars2[i]
        
        xp <- x2[parameter == param, .(val)]
        
        dens <- density(xp$val, adjust = 2, cut = 0)
        dd <- with(dens, data.table(x, y))
        
        sig_figs <- 3
        mx <- median(xp$val)
        if (ci == "hpdi") {
            hpdi <- hdi(xp$val, credMass = 0.95)
            lq <- hpdi[["lower"]]
            uq <- hpdi[["upper"]]
            # list(lq = lq, uq = uq)
        } else {
            lq <- quantile(xp$val, 0.025)
            uq <- quantile(xp$val, 0.975)
        }
        
        if (param == "cov_G_ii") {
            rng <- qnorm(c(0.1, 0.9), 0, sqrt(mx)) |> diff() |> exp()
            rng_str <- paste0(" x", signif(rng, 3))
            print(glue("G cov(inf) = {signif(mx, 3)} ({signif(lq, 3)}, {signif(uq, 3)}), {rng_str}"))
        }
        
        h2_str <- if (startsWith(param, "cov_G")) {
            h2var <- sub("cov_", "h2_", param)
            h2 <- x[, mean(get(h2var))]
            message("h2 in (", signif(x[, quantile(get(h2var), .025)], 3), ", ",
                    signif(x[, quantile(get(h2var), .975)], 3), ")")
            glue(", h<sup>2</sup>={signif(h2, 2)}")
            ""
        } else {
            ""
        }
        
        if (startsWith(param, "cov")) {
            x_min <- 0
            x_max <- 1.0
        } else {
            x_min <- -1
            x_max <- 1
        }
        
        plts[[param]] <- ggplot() +
            geom_line(data = dd, aes(x = x, y = y)) +
            geom_ribbon(data = subset(dd, lq < x & x < uq),
                        aes(x = x, ymin = 0, ymax = y),
                        fill = "red", colour = NA, alpha = 0.5) +
            geom_vline(xintercept = mx, linewidth = 1, colour = "blue") +
            coord_cartesian(xlim = c(x_min, x_max)) +
            labs(x = "value", y = "density",
                 title = glue("{param2}: {signif(mx, 2)} ({signif(lq, 2)}, {signif(uq, 2)}){h2_str}")) +
            theme(plot.title = element_markdown(),
                  text = element_text(size = 10),
                  axis.text.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  # text = element_text(size = 28),
                  legend.position = "none")
    }
    
    plts[["empty"]] = ggplot() + theme_void()
    
    pltG_list <- c("cov_G_ss", "r_G_si", "r_G_sr",
                   "empty", "cov_G_ii", "r_G_ir",
                   "empty", "empty", "cov_G_rr")
    pltE_list <- sub("G", "E", pltG_list)
    
    pltG <- plot_grid(plotlist = plts[pltG_list])
    pltE <- plot_grid(plotlist = plts[pltE_list])
    
    title_plt <- ggplot() + labs(title = title_plot_str)
    
    pltG_cov <- pltG # plot_grid(title_plt, pltG, ncol = 1, rel_heights = c(0.06, 1))
    pltE_cov <- plot_grid(title_plt, pltE, ncol = 1, rel_heights = c(0.06, 1))
    pltGE_cov <- plot_grid(title_plt, plot_grid(pltG, pltE, ncol = 2), ncol = 1, rel_heights = c(0.06, 1))
    
    
    dir_str <- glue("gfx/{data_set}")
    if (!dir.exists(dir_str)) {
        dir.create(dir_str)
    }
    
    pngG_str <- glue("{dir_str}/cov-G-{scenario}.png")
    pngE_str <- glue("{dir_str}/cov-E-{scenario}.png")
    pngGE_str <- glue("{dir_str}/cov-GE-{scenario}.png")
    
    message("Saving to ", pngG_str)
    ggsave(pngG_str, pltG_cov, width = 3000, height = 2000, units = "px")
    
    message("Saving to ", pngE_str)
    ggsave(pngE_str, pltE_cov, width = 3000, height = 2000, units = "px")
    
    message("Saving to ", pngGE_str)
    ggsave(pngGE_str, pltGE_cov, width = 4000, height = 2000, units = "px")
    
    list(G = pltG_cov, E = pltE_cov, GE = pltGE_cov)
}

