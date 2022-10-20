library(data.table)
library(glue)
library(gtools)
library(ggplot2)
library(cowplot)
library(viridis)
library(RColorBrewer)

source("utils.R")

plot_violins <- function(data_set = "fb_1_12_Dx", split_replicates = FALSE) {
    
    # data_set = "fb", split_replicates = FALSE
    
    parameters <- c("beta", "latent_period", "eta_shape", "detection_period",
                    "rho_shape", "recovery_period", "gamma_shape",
                    "donor_i", "donor_l", "donor_d", "donor_r",
                    "trial_i", "trial_l", "trial_d", "trial_r",
                    "sigma")
    
    params2 <- c("beta", "Latent Period (days)", "LP Shape", "Detection Period (days)",
                 "DP Shape", "Recovery Period (days)", "RP Shape",
                 "Donor infectivity", "Donor latency", "Donor detectability", "Donor recoverability",
                 "Trial infectivity", "Trial latency", "Trial detectability", "Trial recoverability",
                 "Group Effect")
    
    gfx_dir <- glue("gfx/{data_set}")
    if (!dir.exists(gfx_dir)) dir.create(gfx_dir)
    
    # don't want shape
    to_keep <- !grepl("_shape", parameters)
    parameters <- parameters[to_keep]
    params2 <- params2[to_keep]
    
    # colours <- viridis(length(parameters))
    # names(colours) <- parameters
    
    # Data wil be stored here in long format
    X <- data.table()
    
    load(list.files(path = glue("results/{data_set}"), full.names = TRUE)[1])
    burn_prop <- if (exists("params")) with(params, burnin / nsample) else 0.2
    burn_prop <- clamp(burn_prop, 0, 1)
    
    # Probably 1:13?
    files <- mixedsort(list.files(path = glue("data/{data_set}"), pattern = "_out"), decreasing = TRUE)
    for (i in seq_along(files)) {
        fi <- sub("_out", "", files[i])
        trace_file <- glue("data/{data_set}/{fi}_out/trace.txt")
        if (file.exists(trace_file)) {
            xi <- fread(trace_file)[(ceiling(burn_prop * .N) + 1):.N]
        } else {
            message(trace_file, " is missing")
            next
        }
        
        # message(fi, " = ", xi[, .N])
        # next
        
        res_file <- glue("results/{data_set}/{fi}.RData")
        if (file.exists(res_file)) {
            load(res_file)
        } else {
            message(glue("{res_file} {fi} is missing"))
            next
        }
        
        # message("Loaded '", fi, "' (", params$label, "), ", params$description)
        message(glue("Loaded '{fi}' ({params$label}), {params$description}"))
        
        for (param in parameters) {
            if (param %in% names(xi)) {
                xi_dt <- data.table(parameter = param,
                                    group = substring(params$label, 1, 2),
                                    label = params$label,
                                    replicate = params$replicate,
                                    value = xi[[param]])
                X <- rbind(X, xi_dt)
            }
        }
    }
    
    # X <- X[group == "s1"]
    
    # for (i in seq_along(parameters)) {
    # X[parameter == parameters[i], parameter := params2[i]]
    # }
    
    X[, `:=`(group = factor(group),
             label = factor(label),
             scenrep = if (split_replicates) {
                 factor(paste0(label, "_", replicate))
             } else {
                 factor(label)
             })]
    
    X_lims <- X[, .(min = min(value),
                    low  = pmin(0, 0.9 * quantile(value, 0.025), 1.1 * quantile(value, 0.025)),
                    high = pmax(0, 0.9 * quantile(value, 0.975), 1.1 * quantile(value, 0.975)),
                    max = max(value)),
                parameter]
    
    setkey(X, parameter)
    setkey(X_lims, parameter)
    
    if (length(unique(X$label)) == 16) {
        levels(X$label) <- list("None" = "s1a", "None" = "s2a",
                                "Lat" = "s1b", "Lat" = "s2b",
                                "Inf" = "s1c", "Inf" = "s2c",
                                "Det" = "s1d", "Det" = "s2d",
                                "Lat+Inf" = "s1e", "Lat+Inf" = "s2e",
                                "Lat+Det" = "s1f", "Lat+Det" = "s2f",
                                "Inf+Det" = "s1g", "Inf+Det" = "s2g",
                                "Lat+Inf+Det" = "s1h", "Lat+Inf+Det" = "s2h")
    } else {
        levels(X$label) <- list("None" = "s1a",
                                "Lat" = "s1b",
                                "Inf" = "s1c",
                                "Det" = "s1d",
                                "Lat+Inf" = "s1e",
                                "Lat+Det" = "s1f",
                                "Inf+Det" = "s1g",
                                "Lat+Inf+Det" = "s1h")
    }
    
    scale_breaks = X[, levels(unique(label))]
    n_breaks = length(scale_breaks)
    
    plts <- list()
    for (i in seq_along(params2)) {
        param <- parameters[i]
        plts[[param]] <-
            ggplot(X[param],
                   aes(x = scenrep, y = value, fill = label)) +
            # geom_violin(scale = "width") +
            geom_boxplot(outlier.shape = NA) +
            geom_hline(yintercept = 0, linetype = "dashed") +
            ylim(X_lims[param, low], X_lims[param, high]) +
            scale_fill_manual("label",
                              breaks = scale_breaks,
                              values = brewer.pal(n_breaks, "Set1")) +
            # values = viridis(n_breaks)) +
            # scale_colour_viridis(discrete = TRUE) +
            labs(x = "Data set", y = params2[i]) +
            scale_x_discrete(limits = rev) +
            coord_flip() +
            theme(text = element_text(size = 16),
                  legend.position = "none")
        
        # png(glue("{gfx_dir}/violin-{param}.png"), width = 2000, height = 1000)
        png(glue("{gfx_dir}/boxplot-{param}.png"), width = 2000, height = 1000)
        print(plts[[param]])
        dev.off()
    }
    
    legend <- get_legend(
        plts$beta +
            guides(color = guide_legend(nrow = 1)) +
            theme(legend.position = "bottom")
    )
    
    plot_periods <- FALSE
    if (plot_periods) {
        png(glue("{gfx_dir}/periods.png"), width = 1000, height = 1000)
        plot_grid(plotlist = plts[c("latent_period", "detection_period", "recovery_period")], nrow = 3)
        dev.off()
    }
    
    plot_donor_FE <- FALSE
    if (plot_donor_FE) {
        png(glue("{gfx_dir}/donor_FEs.png"), width = 1500, height = 1000)
        plot_grid(plotlist = plts[startsWith(names(plts), "donor_")], nrow = 2)
        dev.off()
    }
    
    plot_trial_FE <- FALSE
    if (plot_trial_FE) {
        png(glue("{gfx_dir}/trial_FEs.png"), width = 1500, height = 1000)
        plot_grid(plotlist = plts[startsWith(names(plts), "trial_")], nrow = 2)
        dev.off()
    }
    
    
    plot_all <- FALSE
    if (plot_all) {
        pltg <- plot_grid(plotlist = plts[names(plts) != "sigma"])
        # pdf("foo2.pdf", width = 16, height = 9)
        # png("gfx/fb-violins/all_pars.png", width = 2000, height = 1000)
        if (split_replicates) {
            png(glue("{gfx_dir}/all_pars_split.png"), width = 4000, height = 2000)
        } else {
            png(glue("{gfx_dir}/all_pars.png"), width = 2000, height = 1000)
        }
        print(pltg)
        dev.off()
    } else {
        gc()
        # print(pltg)
        # p1 <- plot_grid(plotlist = plts[c("beta", "latent_period", "detection_period", "donor_i", "donor_l", "donor_d", "trial_l", "trial_d")])
        p1 <- plot_grid(plotlist = plts[names(plts) != "sigma"])
        pltg <- plot_grid(p1, legend, ncol = 1, rel_heights = c(1, 0.1))
        if (split_replicates) {
            png(glue("{gfx_dir}/all_pars_split.png"), width = 4000, height = 2000)
        } else {
            png(glue("{gfx_dir}/all_pars.png"), width = 2000, height = 1000)
        }
        print(pltg)
        dev.off()
    }
    
    file_str <- glue("meta/boxplots-{data_set}.RData") # -{scenario}
    message("saving to ", file_str)
    save(X, plts, pltg, file = file_str)
    
    plts
}


if (FALSE) {
    plts <- plot_violins(data_set = "fb_1_12_Dx_linked", split_replicates = FALSE)
    plts_split <- plot_violins(data_set = "fb_1_12_Dx_linked", split_replicates = TRUE)
}

# plts_split <- plot_violins(data_set = "sim-simple", split_replicates = TRUE)

if (FALSE) {
    x <- plts_split
    legend <- get_legend(
        x$beta +
            guides(color = guide_legend(nrow = 1)) +
            theme(legend.position = "bottom")
    )
    
    plot_set <- c("beta", "latent_period", "detection_period", "recovery_period")
    
    p1 <- plot_grid(plotlist = x[plot_set])
    pltg_split <- plot_grid(p1, legend, ncol = 1, rel_heights = c(1, 0.1))
    pltg_split
}

# save(plts_fb_12_Dx_link, plts_fb_12_Dx_split,
#      file="Thursday1.RData")
