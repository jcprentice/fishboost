{
    library(data.table)
    library(tidyverse)
    library(ggplot2)
    library(ggtext)
    library(cowplot)
    # library(ggh4x)
    source("utils.R")
}

# gg_colour_hue <- function(n) {
#     hues = seq(15, 375, length = n + 1)
#     hcl(h = hues, l = 65, c = 100)[1:n]
# }

model_fit_rmsd <- function(dataset = "fb-test", scens = 0, alt = "", drop_outliers = FALSE) {
    if (FALSE) {
        dataset <- "fb-test"; scens <- c(1,3,4); drop_outliers = FALSE; alt = ""
        dataset <- "sim-base-inf"; scens <- 1:16; drop_outliers = FALSE; alt = ""
    }
    
    message(str_glue("Calculating RMS model fit for '{dataset}'"))
    
    if (drop_outliers) message("- Dropping outliers")
    outliers <- if (drop_outliers) "-drop" else ""
    
    if (str_1st(alt) %notin% c("", "-")) alt <- str_c("-", alt)
    
    km_data <- readRDS(str_glue("datasets/{dataset}/meta/km_data_ps.rds"))
    
    if (all(scens != 0)) {
        # km_data <- km_data[scens]
    }
    
    params <- map(km_data, "params")
    labels <- map_chr(params, "label")
    scenarios <- map_int(params, "scenario")
    
    # Set limits on AUC length
    tmax <- if (str_detect(dataset, "fb")) {
        c(104, 160)
    } else {
        c(200, 200)
    }
    
    # Extract RMS deviance for each id, sire, trial / Tsym, RP
    fit <- map(km_data, \(x) {
        # i <- 1; x <- km_data[[i]]
        x1 <- x$data[, .(id, sire, trial, Tsym, RP, src)] |>
            setorder(id, sire, trial)
        
        if (drop_outliers) {
            drop_ids <- if (km_data[[1]]$data[, max(sire)] == 28) {
                c(3, 8, 22, 25)
            } else {
                c(3, 8, 23, 26)
            } 
            x1 <- x1[sire %notin% drop_ids]
        }
        
        x1[, `:=`(Tsym = as.integer(fifelse(Tsym > tmax[trial], tmax[trial], ceiling(Tsym), tmax[trial])),
                  RP   = as.integer(fifelse(RP   > tmax[trial], tmax[trial], ceiling(RP),   tmax[trial])))]
        
        # x1[, `:=`(Tsym = fifelse(Tsym > tmax[trial], tmax[trial], Tsym, tmax[trial]),
        #           RP   = fifelse(RP   > tmax[trial], tmax[trial], RP,   tmax[trial]))]
        
        x2 <- x1[, .(Tsym = c(0, sort(Tsym, na.last = TRUE)),
                     RP = c(0, sort(RP, na.last = TRUE)),
                     survival = seq(1, 0, length = .N + 1),
                     src = src[[1]]),
                 .(id, sire, trial)] |>
            melt(measure.vars = c("Tsym", "RP"))
        
        x3 <- x2[, {
            # message(str_glue("id == {id[[1]]} & sire == {sire[[1]]} & trial == {trial[[1]]} & variable == \"{variable[[1]]}\""))
            times <- seq(0, tmax[trial])
            if (last(value) < tmax[trial]) {
                value <- c(value, tmax[trial])
                survival <- c(survival, last(survival))
            }
            f <- approxfun(value, survival,
                           method = "constant",
                           ties = list("ordered", max))
            s1 <- f(times)
            s1[is.na(s1)] <- 0
            list(time = times, survival = s1, src = src[[1]])
        }, .(id, sire, trial, variable)]
        
        x4 <- x3[, dev := survival - survival[[.N]], time][id != id[[.N]]]
        
        x4[, .(all  = sqrt(mean(dev^2)),
               Tsym = sqrt(mean(dev[variable == "Tsym"]^2)),
               RP   = sqrt(mean(dev[variable == "RP"]^2))),
               # mad_all   = mean(abs(dev)),
               # mad_Tsym  = mean(abs(dev[variable == "Tsym"])),
               # mad_RP    = mean(abs(dev[variable == "RP"]))),
           .(id, sire, trial)] |>
            melt(measure.vars = c("all", "Tsym", "RP"))
    }) |>
        rbindlist(idcol = "scen")
    
    fit[, scen := ordered(labels[scen], levels = str_sort(labels, numeric = TRUE))]
    
    # fit[, `:=`(scen = factor(str_c("s", scen), levels = labels),
    #            id = format(id),
    #            sire = format(sire),
    #            trial = format(trial))]
    setorder(fit, scen, id, sire, trial, variable)
    
    fit_wide <- fit[, .(value = mean(value)), .(scen, variable)] |>
        dcast(scen ~ variable)
    fit_wide
    
    
    fit_plt <- ggplot(fit,
                      aes(x = scen, y = value, colour = variable)) +
        stat_summary(geom = "errorbar", fun.data = mean_se,
                     position = position_dodge(0.1),
                     width = 0.2) +
        stat_summary(geom = "point", fun = mean,
                     position = position_dodge(0.1)) +
        # geom_point(show.legend = FALSE,
        #            position = position_dodge(0.25)) +
        scale_y_continuous(limits = ~ range(.x, 0)) +
        scale_colour_discrete("Measurements",
                            breaks = c("all", "Tsym", "RP"),
                            # values = gg_colour_hue(3),
                            labels = c("Combined", "Tsym", "RP")) +
        labs(x = "Scenario",
             y = "RMS of deviation",
             title = "RMS of deviation between data and model output",
             subtitle = str_glue("Dataset: '{dataset}'")) +
        theme_classic() +
        theme(plot.title = element_text(size = 16),
              plot.background = element_rect(fill = "white"))
    fit_plt
    
    mf_dir <- str_glue("datasets/{dataset}/gfx/model_fit")
    if (!dir.exists(mf_dir)) {
        message("- mkdir ", mf_dir)
        dir.create(mf_dir)
    }
    mf_base <- str_glue("{mf_dir}/{dataset}-model_fit-rmv_dev")
    
    plt_str <- str_glue("{mf_base}{outliers}{alt}.png")
    ggsave(plt_str, fit_plt, width = 8, height = 5)
    
    families_plt <- ggplot(fit, aes(x = sire, y = value, colour = variable)) +
        stat_summary(geom = "errorbar", fun.data = mean_se,
                     position = position_dodge(0.1),
                     width = 0.2) +
        stat_summary(geom = "point", fun = mean,
                     position = position_dodge(0.1)) +
        # geom_point(show.legend = FALSE,
        #            position = position_dodge(0.25)) +
        scale_colour_discrete("Values",
                            breaks = c("all", "Tsym", "RP"),
                            # values = gg_colour_hue(3),
                            labels = c("Combined", "Tsym", "RP")) +
        labs(x = "Family",
             y = "RMS of deviation",
             title = "RMS of deviation between data and model output, per family",
             subtitle = str_glue("Dataset: '{dataset}'")) +
        facet_wrap(~ scen, ncol = 4) +
        theme_bw() +
        theme(legend.position = "bottom")
    families_plt
    
    plt_str <- str_glue("{mf_base}_by_family{outliers}{alt}.png")
    ggsave(plt_str, families_plt, width = 12, height = 8)
    
    list(fit = fit_wide,
         fit_plt = fit_plt,
         families_plt = families_plt)
}

out_fb_final <- model_fit_rmsd("fb-final")
out_fb_final <- model_fit_rmsd("fb-final2")
# out_fb_final_drop <- model_fit_rmsd("fb-final", drop_outliers = TRUE)
out_fb_lp <- model_fit_rmsd("fb-lp")
out_fb_lp2 <- model_fit_rmsd("fb-lp2")
# out_fb_lp_drop <- model_fit_rmsd("fb-lp", drop_outliers = TRUE)
out_fb_simple <- model_fit_rmsd("fb-simple")
out_fb_donors <- model_fit_rmsd("fb-donors")
# out_fb_donors_drop <- model_fit_rmsd("fb-donors", drop_outliers = TRUE)
out_sbi <- model_fit_rmsd("sim-base-inf", 0, "", FALSE)
out_fb_test <- model_fit_rmsd("fb-test", 0, "")
