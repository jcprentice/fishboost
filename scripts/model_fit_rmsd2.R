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

model_fit_rmsd2 <- function(dataset = "fb-final", scens = 1:8,
                            alt = "", opts = NULL) {
    # dataset <- "fb-final"; drop_outliers = FALSE; alt = ""
    # dataset <- "sim-base-inf"; scens <- 1:14; drop_outliers = FALSE; alt = ""
    message(str_glue("Calculating RMS model fit for '{dataset}'"))
    
    if ("drop_outliers" %in% opts) message(" - Dropping outliers")
    outliers <- if ("drop_outliers" %in% opts) "-drop" else ""
    
    if (str_1st(alt) %notin% c("", "-")) alt <- str_c("-", alt)
    
    km_data <- readRDS(str_glue("datasets/{dataset}/meta/km_data_ps.rds"))
    
    if (all(scens != 0)) {
        km_data <- km_data[scens]
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
        x1 <- x$data[!is.na(sire), .(id, sire, trial, Tsym, RP, src)] |>
            setorder(id, sire, trial)
        
        # Round up times to represent polling
        x1[, `:=`(Tsym = as.integer(fifelse(Tsym > tmax[trial], tmax[trial], ceiling(Tsym), tmax[trial])),
                  RP   = as.integer(fifelse(RP   > tmax[trial], tmax[trial], ceiling(RP),   tmax[trial])))]
        
        # x1[, `:=`(Tsym = fifelse(Tsym > tmax[trial], tmax[trial], Tsym, tmax[trial]),
        #           RP   = fifelse(RP   > tmax[trial], tmax[trial], RP,   tmax[trial]))]
        
        if ("drop_outliers" %in% opts) {
            drop_ids <- if (max(x1$sire, na.rm = TRUE) == 28) {
                c(3, 8, 22, 25)
            } else {
                c(3, 8, 23, 26)
            }
            x1 <- x1[sire %notin% drop_ids]
        }
        
        # Create survival curves by adding a start (0, 1), followed by times and
        # survival (proportion remaining), then melt
        x2 <- x1[, .(Tsym = c(0, sort(Tsym, na.last = TRUE)),
                     RP = c(0, sort(RP, na.last = TRUE)),
                     survival = seq(1, 0, length = .N + 1),
                     src = src[[1]]),
                 .(id, sire, trial)] |>
            melt(measure.vars = c("Tsym", "RP"))
        
        # Use `approxfun()` to resample the points as a step function with all
        # times included
        x3 <- x2[, {
            # message(str_glue("id == {id[[1]]} & sire == {sire[[1]]} & trial == {trial[[1]]} & variable == \"{variable[[1]]}\""))
            times <- seq(0, tmax[trial])
            if (last(value) < tmax[trial]) {
                value <- c(value, tmax[trial])
                survival <- c(survival, last(survival))
            }
            f <- approxfun(value, survival,
                           method = "constant",
                           ties = list("ordered", min))
            s1 <- f(times)
            s1[is.na(s1)] <- 0
            list(time = times, survival = s1, src = src[[1]])
        }, .(id, sire, trial, variable)]
        
        
        # Calculate the mean survival curve across all sims
        x4 <- x3[, .(survival = mean(survival)), .(sire, trial, variable, time, src)] |>
            dcast(... ~ src, value.var = "survival") |>
            setorder(variable, trial, sire, time) |>
            setkey(NULL)
        x4[, `:=`(dev = fb - sim, fb = NULL, sim = NULL)]
        
        rms <- function(x) x^2 |> mean() |> sqrt()
        
        x5 <- x4[, .(all  = rms(dev),
                     Tsym = rms(dev[variable == "Tsym"]),
                     RP   = rms(dev[variable == "RP"])),
                 .(sire, trial)]
        x5
    }) |>
        rbindlist(idcol = "scen")
    
    fit[, scen := ordered(labels[scen], levels = str_sort(labels, numeric = TRUE))]
    
    mfit <- melt(fit, measure.vars = c("all", "Tsym", "RP"))
    
    fit_plt <- ggplot(mfit,
                      aes(x = scen, y = value, colour = variable)) +
        stat_summary(geom = "errorbar",
                     fun.data = mean_se,
                     fun.args = list(mult = 1.96),
                     position = position_dodge(0.1),
                     width = 0.2) +
        stat_summary(geom = "point",
                     fun = mean,
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
             title = "RMS of deviation of mean model curves vs data",
             subtitle = str_glue("Dataset: '{dataset}'")) +
        theme_classic() +
        theme(plot.title = element_text(size = 16),
              plot.background = element_rect(fill = "white"))
    fit_plt
    
    mf_dir <- str_glue("datasets/{dataset}/gfx/model_fit")
    if (!dir.exists(mf_dir)) dir.create(mf_dir)
    mf_base <- str_glue("{mf_dir}/{dataset}-model_fit-rmv_dev2")
    
    plt_str <- str_glue("{mf_base}{outliers}{alt}.png")
    ggsave(plt_str, fit_plt, width = 8, height = 5)
    
    list(fit = fit,
         fit_plt = fit_plt)
}

