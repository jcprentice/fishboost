{
    library(data.table)
    library(purrr)
    library(stringr)
    library(ggplot2)
    library(ggtext)
    library(cowplot)
    # library(ggh4x)
}

model_fit_auc <- function(dataset = "fb-final", drop_outliers = FALSE) {
    # dataset <- "fb-final"; drop_outliers = FALSE
    
    message(str_glue("Calculating AUC model fit for '{dataset}'"))
    
    if (drop_outliers) message("- Dropping outliers")
    outliers <- if (drop_outliers) "-drop" else ""
    
    km_data <- readRDS(str_glue("datasets/{dataset}/meta/km_data_ps.rds"))
    
    params <- map(km_data, "params")
    labels <- map_chr(params, "label")
    
    # Set limits on AUC length
    ttest <- 1
    tmax <- c(104, 160) * ttest
    
    # Extract AUC for each id, sire, trial / Tsym, RP
    AUC <- map(km_data, \(x) {
        # x <- km_data[[1]]
        x1 <- x$data[!is.na(sire), .(id, sire = as.factor(sire),
                                     trial = as.factor(trial), Tsym, Tdeath)] |>
            setorder(id, sire, trial)
        
        if (drop_outliers) {
            drop_ids <- if (km_data[[1]]$data[, max(sire)] == 28) {
                c(3, 8, 22, 25)
            } else {
                c(3, 8, 23, 26)
            } 
            x1 <- x1[sire %notin% drop_ids]
        }
        
        # TODO: Sensitivity Analysis on the cutoff point
        if (!is.na(ttest)) {
            x1[, `:=`(Tsym   = fifelse(Tsym   < tmax[trial], Tsym,   tmax[trial]),
                      Tdeath = fifelse(Tdeath < tmax[trial], Tdeath, tmax[trial]))]
        }
        x1[, `:=`(RP = Tdeath - Tsym, Tdeath = NULL)]
        x1[, `:=`(Tsym = sort(Tsym, na.last = TRUE),
                  RP = sort(RP, na.last = TRUE),
                  survival = seq(1, 0, length = .N)),
           .(id, sire, trial)]
        
        # Note: mean(x) != sum(x) / N, when x contains NAs. We want to get the
        # AUC from a curve that starts at (0,1) and steps down with each event.
        x1[, .(Tsym = sum(diff(c(0, Tsym)) * survival, na.rm = TRUE),
               RP = sum(diff(c(0, RP)) * survival, na.rm = TRUE)),
           .(id, sire, trial)]
    })
    
    # Data needs to be in long format for plotting, the X value is the single
    # experiment value, the Y value is each simulation value.
    AUC_long <- map(AUC, \(x) {
        # x <- AUC[[1]]
        x1 <- melt(x, measure.vars = c("Tsym", "RP"), value.name = "sim")
        
        x2 <- x1[, obs := sim[[.N]], .(sire, trial, variable)][id != id[[.N]]]
        
        # x2[, .(sim = mean(sim), obs = mean(obs)), .(sire, trial, variable)]
        x2
    })
    
    # Get the values of the model fits, and rank them
    {
        AUC_long_all <- AUC_long |>
            rbindlist(idcol = "scenario")
        AUC_long_all[, `:=`(scenario = str_c("s", scenario),
                            scen_var = str_c("s", scenario, "_", variable))]
        
        fit <- AUC_long_all[, {
            tmp1 <- cor.test(sim, obs)
            tmp2 <- cor.test(sim[variable == "Tsym"], obs[variable == "Tsym"])
            tmp3 <- cor.test(sim[variable == "RP"],   obs[variable == "RP"])
            list(all_est  = tmp1$estimate, all_low  = tmp1$conf.int[[1]], all_high  = tmp1$conf.int[[2]],
                 Tsym_est = tmp2$estimate, Tsym_low = tmp2$conf.int[[1]], Tsym_high = tmp2$conf.int[[2]],
                 RP_est   = tmp3$estimate, RP_low   = tmp3$conf.int[[1]], RP_high   = tmp3$conf.int[[2]])},
            scenario]
        # fit <- AUC_long_all[, .(cor = cor(sim, obs),
        #                         cor_Tsym = cor(sim[variable == "Tsym"], obs[variable == "Tsym"]),
        #                         cor_RP = cor(sim[variable == "RP"], obs[variable == "RP"])),
        #                     scenario]
        (fit[, map(.SD, \(x) if (is.numeric(x)) round(100*x, 1) else x)]
            [, .(scenario,
                 all  = str_c(all_est,  " (", all_low,  ", ", all_high,  ")"),
                 Tsym = str_c(Tsym_est, " (", Tsym_low, ", ", Tsym_high, ")"),
                 RP   = str_c(RP_est,   " (", RP_low,   ", ", RP_high,   ")"))])
    }
    
    
    # fit_long <- fit |> melt(id.vars = "scenario")
    fit_long <- fit |> melt(id.vars = "scenario",
                            measure.vars = measure(variable, range, sep = "_")) |>
        dcast(scenario + variable ~ range)
    fit_long[, `:=`(scenario = ordered(scenario, levels = str_sort(unique(scenario), numeric = TRUE)),
                    variable = factor(variable, levels = c("all", "Tsym", "RP")))]
    setorder(fit_long, scenario, variable)
    
    fit_plt <- ggplot(fit_long,
                      aes(x = scenario, colour = variable)) +
        geom_hline(yintercept = 0, linetype = "dashed", colour = "black") +
        # geom_segment(aes(y = low, yend = high), position = position_dodge(0.1)) +
        geom_errorbar(aes(ymin = low, ymax = high),
                      position = position_dodge(0.2),
                      width = 0.2) +
        geom_point(aes(y = est),
                   position = position_dodge(0.2)) +
        expand_limits(y = c(-0.3, 1)) +
        scale_colour_discrete("Measurements",
                              breaks = c("all", "Tsym", "RP"),
                              labels = c("Combined", "Tsym", "RP")) +
        labs(x = "Scenario",
             y = "Correlation",
             title = "Correlation (with 95% CI) between true and simulated family AUCs",
             subtitle = str_glue("Dataset: '{dataset}'")) +
        theme_classic() +
        theme(plot.title = element_text(size = 16))
    fit_plt
    
    mf_dir <- str_glue("datasets/{dataset}/gfx/model_fit")
    if (!dir.exists(mf_dir)) {
        message("- mkdir ", mf_dir)
        dir.create(mf_dir)
    }
    
    ggsave(str_glue("{mf_dir}/{dataset}-model_fit-auc{outliers}.png"),
           fit_plt, width = 8, height = 5)
    
    
    
    AUC_plts <- map2(AUC_long, labels, \(x, label) {
        # x <- AUC_long[[1]]; label <- labels[[1]]
        ggplot(x, aes(obs, sim, colour = sire)) +
            stat_summary(geom = "errorbar", fun.data = mean_se,
                         show.legend = FALSE) +
            stat_summary(geom = "point", fun = mean,
                         show.legend = FALSE) +
            geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
            geom_smooth(data = x, aes(obs, sim),
                        inherit.aes = FALSE,
                        method = lm, se = FALSE) +
            labs(x = "AUC (data)",
                 y = "AUC (simulation)",
                 title = str_glue("Scenario {label}")) +
            expand_limits(x = 0, y = 0) +
            facet_wrap(~ variable,
                       scales = "free",
                       labeller = labeller(variable = c(
                           "Tsym" = "Time to first symptoms",
                           "RP" = "Time from symptoms to death"))) +
            theme_bw()
    })
    
    title_plt <- ggplot() + labs(title = str_glue("Dataset: '{dataset}'")) + theme_classic()
    
    
    AUC_plts_grid <- plot_grid(title_plt,
                               plot_grid(plotlist = AUC_plts),
                               ncol = 1, rel_heights = c(0.05, 1))
    AUC_plts_grid
    
    
    ggsave(str_glue("{mf_dir}/{dataset}-model_fit-auc_by_family{outliers}.png"),
           AUC_plts_grid, width = 18, height = 12)
    
    list(fit = fit,
         fit_plt = fit_plt,
         AUC_plts = AUC_plts_grid)
}

out_fb_final <- model_fit_auc("fb-final")
out_fb_final2 <- model_fit_auc("fb-final2")
# out_fb_final_drop <- model_fit_auc("fb-final", drop_outliers = TRUE)
out_fb_lp <- model_fit_auc("fb-lp")
out_fb_lp2 <- model_fit_auc("fb-lp2")
out_fb_simple <- model_fit_auc("fb-simple")
out_fb_donors <- model_fit_auc("fb-donors")
