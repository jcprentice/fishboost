library(data.table)
library(ggplot2)
library(cowplot)


plot_km <- function() {
    x <- readRDS("fb_data/fb_data12.rds")
    x <- x[sdp == "progeny"]
    
    # censor final times
    # x[, Trec := fifelse(Trec == max(Trec, na.rm = TRUE), NA_real_, Trec), trial]
    
    
    y <- melt(x[, .(Trial = trial, group = sire, Tsym, Trec = Trec - Tsym)],
              measure.vars = c("Tsym", "Trec"),
              variable.name = "var",
              value.name = "time")
    
    y <- y[order(Trial, group, var, time)]
    
    y[, `:=`(Trial = factor(Trial),
             group = factor(group),
             var = factor(var))]
    
    y1 <- y[, .(time = c(0, sort(time, na.last = TRUE))),
            by = .(Trial, group, var)]
    y1[var == "Tsym", survival := 1 - seq(0, 1, length.out = .N),
       by = .(Trial, group)]
    y1[var == "Trec" & !is.na(time), `:=`(survival = 1 - seq(0, 1, length.out = .N), N = .N),
       by = .(Trial, group)]
    
    plt <- list()
    
    plt[["Tsym"]] <- ggplot(data = y1[var == "Tsym"],
                            aes(x = time,
                                y = survival,
                                group = group,
                                colour = Trial)) +
        # geom_line(size = 2, alpha = 0.5) +
        geom_step(size = 1.5, alpha = 0.5) +
        labs(title = "Time from start to first symptoms", x = "Time (days)", y = "Proportion") +
        theme(legend.position = "none",
              text = element_text(size = 48))

    plt[["Trec"]] <- ggplot(y1[var == "Trec"],
                            aes(x = time,
                                y = survival,
                                group = group,
                                colour = Trial)) +
        # geom_line(size = 2, alpha = 0.5) +
        geom_step(size = 1.5, alpha = 0.5) +
        labs(title = "Time from symptoms to death", x = "Time (days)", y = "Proportion") +
        theme(legend.position = c(0.8, 0.85),
              legend.key.width = unit(0.3, "npc"),
              text = element_text(size = 48))
    
    plts <- plot_grid(plotlist = plt)

    ggsave("gfx/km_plots.png", plts, width = 2000, height = 1000, units = "px")

    plts
}

plot_km()
