source("add_latency.R")

pop <- data.table(id = 1:20, sdp = "progeny", group = rep(1:4, each = 5))
pop[, Tinf := runif(.N, 0, 10)]
pop[pop[, .(I = .I[1]), group]$I, Tinf := 0]
pop[, Tsym := Tinf + runif(.N, 0, 1)]
pop[, Trec := Tsym + runif(.N, 0, 1)]
pop[Tinf != 0 & sample(c(T,F), .N, replace = TRUE, prob = c(0.1, 0.9)), c("Tinf", "Tsym", "Trec") := NA]
pop[]


library(gridExtra)

{
    lp <- add_latency1(pop)
    lp

    pop[, Tinf2 := pmax(Tsym - lp, 0)]


    pop1 <- pop[!is.na(Tsym) & group %in% 2:4, .(group, Tinf, Tinf2, Tsym, Trec)]
    setkey(pop1, group, Tinf, Trec)
    pop1[, id := .N - .I + 1]
}


{
    plt1 <- ggplot(pop1) +
        geom_segment(aes(x = Tinf, xend = Tsym, y = id, yend = id, colour = "Exposed")) +
        geom_segment(aes(x = Tsym, xend = Trec, y = id, yend = id, colour = "Infectious")) +
        geom_point(aes(x = Tinf, y = id, colour = "Tinf")) +
        geom_point(aes(x = Tsym, y = id, colour = "Tsym")) +
        geom_point(aes(x = Trec, y = id, colour = "Trec")) +
        labs(title = "SEIR with Tinfs", x = "Time (days)", y = "id") +
        scale_colour_manual("Compartments",
                            breaks = c("Exposed", "Infectious", "Tinf", "Tsym", "Trec"),
                            values = c("blue", "red", "blue", "purple", "red"),
                            guide = guide_legend(override.aes = list(
                                linetype = c(rep("solid", 2), rep("blank", 3)),
                                shape = c(rep(NA, 2), rep(16, 3))
                            ))) +
        theme(legend.position = c(0.85, 0.75))

    plt2 <- ggplot(pop1) +
        geom_segment(aes(x = Tsym, xend = Trec, y = id, yend = id, colour = "Infectious")) +
        geom_point(aes(x = Tsym, y = id, color = "Tsym")) +
        geom_point(aes(x = Trec, y = id, color = "Trec")) +
        labs(title = "SEIR with missing Tinfs", x = "Time (days)", y = "id") +
        scale_colour_manual("Compartments",
                            breaks = c("Infectious", "Tsym", "Trec"),
                            values = c("red", "purple", "red"),
                            guide = guide_legend(override.aes = list(
                                linetype = c("solid", rep("blank", 2)),
                                shape = c(NA, rep(16, 2))
                            ))) +
        theme(legend.position = c(0.85, 0.8))

    plt <- grid.arrange(plt1, plt2, ncol = 2)
}

library(ggplot2)
library(gridExtra)
library(cowplot)
library(viridis)

source("import_fishboost_data.R")
fb_end <- fb_data[, max(Trec, na.rm = TRUE)]


pdf("FB_plots.pdf", paper = "a4")
for (i in 0:11) {
    gps <- c(1+6*i, 6+6*i)
    fb1 <- fb_data[group %in% seq.int(gps[1], gps[2])]
    fb1[, group := sprintf("%02d", group)]


    fb1[order(-group, -Tsym), id := .I]
    setkey(fb1, id)
    fb1


    plt <- ggplot(fb1, aes(colour = group)) +
        coord_cartesian(xlim = c(0, fb_end)) +
        # geom_segment(aes(x = Tinf, xend = Tsym, y = id, yend = id, colour = "Exposed")) +
        geom_segment(aes(x = Tsym, xend = Trec, y = id, yend = id)) +
        # geom_segment(aes(x = Trec, xend = fb_end, y = id, yend = id, colour = "Removed")) +
        # geom_point(aes(x = Tinf, y = id, colour = "Tinf")) +
        geom_point(aes(x = Tsym, y = id, colour = "t(I)")) +
        geom_point(aes(x = Trec, y = id, colour = "t(R)")) +
        labs(title = paste0("Fishboost: groups ", gps[1], "-", gps[2]),
             x = "Time (days)", y = "Individual") +
        scale_colour_manual("Groups",
                            breaks = waiver(),
                            values = c(viridis(6), "blue", "red"),
                            guide = guide_legend(override.aes = list(
                                linetype = c(rep("solid", 6), rep("blank", 2)),
                                shape = c(rep(NA, 6), rep(16, 2))))) +
        theme(aspect.ratio = sqrt(2))

    print(plt)
}
dev.off()


