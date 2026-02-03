{
    library(data.table)
    library(ggplot2)
    library(purrr)
    library(stringr)
    
    source("fes_to_vals.R")
}

draw_event_times <- function(dataset = "fb-test", scen = 1, rep = 1) {
    # dataset <- "testing"; scen <- 1; rep <- 1
    
    f <- str_glue("datasets/{dataset}/results/scen-{scen}-{rep}.rds")
    res <- readRDS(f)
    
    fb_times <- (
        res$popn
        [sdp == "progeny"]
        [, `:=`(trial = str_c("Tr", trial),
                donor = fifelse(donor == 1, "Don", "Rec"))]
        [, .(# Tinf = mean(Tinf, na.rm = TRUE),
             # Tinc = mean(Tinc, na.rm = TRUE),
             Tsym = mean(Tsym, na.rm = TRUE),
             Tdeath = mean(Tdeath, na.rm = TRUE),
             .N),
            .(trial, donor)]
    )
    
    
    pe <- res$parameter_estimates[str_detect(parameter, "period"),
        .(parameter, mean)] # , hdi95min, hdi95max
    
    pe[, parameter := str_replace_all(
        parameter,
        c("latent_period" = "LP",
          "detection_period" = "DP",
          "removal_period" = "RP",
          "," = "_")
    )]
    pe[, `:=`(period = str_split_i(parameter, "_", 1),
              trial = str_split_i(parameter, "_", 2),
              donor = str_split_i(parameter, "_", 3))]
    pe[, parameter := NULL]
    setcolorder(pe, c("trial", "donor", "period", "mean"))
    setorder(pe, trial, donor)
    
    
    pe[, `:=`(start = cumsum(data.table::shift(mean, 1, 0)),
              end = cumsum(mean))]
    
    x <- merge(pe, fb_times, by = c("trial", "donor"))
    
    x[, `:=`(start = start + Tsym - start[period == "RP"],
             end = end + Tsym - start[period == "RP"]),
      .(trial, donor)]
    x[, c("mean", "Tsym") := NULL]
    
    x[, name := forcats::fct_rev(ordered(str_c(trial, "_", donor)))]
    x[, group := as.integer(name)]
    period_lvls <- c("LP", "DP", "RP")
    x[, period := ordered(period, levels = period_lvls)]
    
    setcolorder(x, c("name", "group"))
    
    plt <- ggplot(x, aes(y = name)) +
        geom_rect(aes(xmin = start, xmax = end,
                      ymin = group - 0.2, ymax = group + 0.2,
                      fill = period)) +
        geom_vline(xintercept = 0,
                   linetype = "dashed", linewidth = 1) +
        labs(x = "Time (days)",
             y = "Group",
             fill = "Period",
             title = str_glue("{dataset} / {scen}-{rep}: {res$params$description}")) +
        theme_bw() +
        theme(plot.title = element_text(size = 8))
    
    gfx_dir <- str_glue("datasets/{dataset}/gfx/events")
    if (!dir.exists(gfx_dir)) {
        message("- mkdir ", gfx_dir)
        dir.create(gfx_dir)
    }
    
    fp <- str_glue("{gfx_dir}/{dataset}-{scen}-{rep}-events.png")
    ggsave(fp, plt, width = 6, height = 3)
    
    plt
}

dataset <- "sim-test"; scen <- 1;

reps <- list.files(str_glue("datasets/{dataset}/results"),
                   str_glue("scen-{scen}")) |>
    str_remove_all("scen-|\\.rds") |>
    str_split_i("-", 2) |>
    as.integer()

walk(reps, ~ draw_event_times(dataset, scen, .x))
