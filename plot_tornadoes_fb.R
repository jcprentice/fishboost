library(data.table)
library(ggplot2)
library(cowplot)
library(extrafont)
library(coda)
library(glue)

plot_tornadoes_fb <- function(scen = 1, data_set = "fb") {

    data_dir <- glue("data/{data_set}")
    results_dir <- sub("data", "results", data_dir)

    data_file <- glue("{results_dir}/scen-fb-{scen}-1.RData")
    if (!file.exists(data_file)) {
        message("file doesen't exist")
        return(NULL)
    }
    load(data_file)
    description <- params$description

    tracetxt <- glue("{data_dir}/scen-fb-{scen}-1_out/trace.txt")
    message(glue("Loading '{tracetxt}'"))
    x <- fread(tracetxt)

    parameters <- names(x)
    parameters <- parameters[!startsWith(parameters, "Group effect") & !startsWith(parameters, "L_")]
    parameters <- parameters[!parameters %in% c("state", "Prior", "Number infected")]

    params2 <- parameters
    params2 <- sub("latent_period", "Latent Period (days)", params2)
    params2 <- sub("eta_shape", "LP shape", params2)
    params2 <- sub("detection_period", "Detection Period (days)", params2)
    params2 <- sub("rho_shape", "DP shape", params2)
    params2 <- sub("recovery_period", "Recovery Period (days)", params2)
    params2 <- sub("gamma_shape", "RP shape", params2)
    params2 <- sub("sigma", "Group Effect", params2)
    params2 <- sub("donor", "Donor FE", params2)
    params2 <- sub("trial", "Trial FE", params2)
    params2 <- sub("_l", " (latency)", params2)
    params2 <- sub("_i", " (infectivity)", params2)
    params2 <- sub("_d", " (detectability)", params2)
    params2 <- sub("_r", " (recoverability)", params2)

    x <- x[1e3:.N, ..parameters]

    xlims <- data.table(parameter = parameters,
                        left  = unlist(x[, lapply(.SD, min)]),
                        right = unlist(x[, lapply(.SD, max)]),
                        key = "parameter")
    xlims["beta",             `:=`(left = 0, right = 0.25)]
    xlims["detection_period", `:=`(left = 0, right = 50)]
    xlims["latent_period",    `:=`(left = 0, right = 100)]
    xlims["eta_shape",        `:=`(left = 0, right = 3)]
    xlims["gamma_shape",      `:=`(left = 0, right = 1.8)]

    # title_str <- paste0("Fitting ", toupper(model), " model")
    title_plt <- ggplot() +
        labs(title = description) +
        theme(text = element_text(size = 24))

    plots <- list()
    for (i in seq_along(parameters)) {
        param <- parameters[i]
        p2 <- params2[i]
        mu_param <- mean(x[[param]])
        plots[[param]] <- ggplot(x, aes_string(x = param)) +
            geom_histogram(aes(y = ..density..),
                           fill = "red",
                           bins = 100,
                           position = "identity",
                           alpha = 0.5) +
            labs(x = p2, y = "Density") +
            coord_cartesian(xlim = c(xlims[param, left], xlims[param, right])) +
            theme(text = element_text(size = 24))
            # geom_text(x = mu_param, y = 30, label = signif(mu_param, 3))

    }

    plt_pars <- plot_grid(
        title_plt,
        plot_grid(plotlist = plots),
        ncol = 1, rel_heights = c(0.06, 1))

    ggsave(glue("gfx/fb-scens/{scen}-pars.png"),
           plt_pars,
           width = 2000, height = 1000, units = "px")

    plt_pars
}


tornado_plots_fb_trials <- function(scen = 1) {

    data_dir <- paste0("data/fb")
    results_dir <- sub("data", "results", data_dir)

    load(glue("{results_dir}/scen-fb-{scen}-1.RData"))
    description <- params$description

    tracetxt1 <- glue("{data_dir}/scen-fb-{scen + 0}-1_out/trace.txt")
    tracetxt2 <- glue("{data_dir}/scen-fb-{scen + 1}-1_out/trace.txt")
    tracetxt3 <- glue("{data_dir}/scen-fb-{scen + 2}-1_out/trace.txt")
    message(glue("Loading '{tracetxt1}'"))
    x1 <- fread(tracetxt1)
    x2 <- fread(tracetxt2)
    x3 <- fread(tracetxt3)

    parameters <- names(x1)
    parameters <- parameters[!startsWith(parameters, "Group effect") & !startsWith(parameters, "L_")]
    parameters <- parameters[!parameters %in% c("state", "Prior", "Number infected")]


    x1 <- x1[, ..parameters]
    x2 <- x2[, ..parameters]
    x3 <- x3[, ..parameters]

    trial_levels <- c("Both", "1 only", "2 only")
    x1[, trial := factor("Both", trial_levels)]
    x2[, trial := factor("1 only", trial_levels)]
    x3[, trial := factor("2 only", trial_levels)]

    x <- rbind(x1, x2, x3)

    # title_str <- paste0("Fitting ", toupper(model), " model")
    title_plt <- ggplot() + labs(title = description)

    plots <- list()
    for (param in parameters) {
        plots[[param]] <- ggplot(x, aes_string(x = param)) +
            geom_histogram(aes(y = ..density.., fill = trial),
                           position = "identity",
                           bins = 100, alpha = 0.5) +
            theme(text = element_text(size = 12))
    }

    plt_pars <- plot_grid(
        title_plt,
        plot_grid(plotlist = plots, nrow = 3),
        ncol = 1, rel_heights = c(0.06, 1))

    ggsave(glue("gfx/fb-trials{scen}-pars.png"),
           plt_pars,
           width = 2000, height = 1000, units = "px")

    plt_pars
}

plt <- plot_tornadoes_fb(scen = 1)
print(plt)
# tornado_plots_fb_trials(scen = 1)
