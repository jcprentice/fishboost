library(data.table)
library(ggplot2)
library(ggeasy)


plot_km <- function() {
    x <- readRDS("fb_data/fb_12.rds")[
        sdp == "progeny", .(trial = factor(trial),
                            # this ensures groups are properly separated
                            group = factor(sprintf("%02d_%02d", trial, sire)),
                            Tsym,
                            RP = Tdeath - Tsym)] |>
        melt(measure.vars = c("Tsym", "RP"),
             variable.name = "var",
             value.name = "time") |>
        setorder(group)
    
    tmax <- c(104, 160)
    x[, time := fifelse(time >= tmax[trial], 161, time, NA)]
    
    y <- x[, .(time = c(0, sort(time, na.last = TRUE)), N = .N),
           by = .(trial, group, var)]
    
    y[, survival := seq(1, 0, length.out = .N),
      by = .(var, trial, group)]
    
    
    plt <- ggplot(data = y,
                  aes(x = time,
                      y = survival,
                      group = group,
                      colour = trial)) +
        geom_step() +
        scale_colour_manual("Trial",
                            breaks = 1:2,
                            labels = c("Data Trial 1", "Data Trial 2"),
                            values = c("blue", "red")) +
        coord_cartesian(xlim = c(0, 160)) +
        # lims(x = c(0, 160), y = 0:1) +
        labs(x = "Time (days)",
             y = "Proportion",
             title = NULL) +
        facet_wrap(. ~ var,
                   nrow = 1,
                   labeller = labeller(
                       var = c(Tsym = "Proportion of family with no symptoms vs time",
                               RP   = "Proportion of family surviving vs time"))) +
        theme_bw() +
        easy_all_text_size(16) +
        theme(panel.background = element_blank(),
              legend.position = "bottom")
    plt
    
    ggsave("gfx/km_plots2.png",
           plt, width = 12, height = 6)

    plt
}
