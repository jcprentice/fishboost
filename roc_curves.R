library(data.table)
library(pROC)
library(ggplot2)
library(cowplot)
# library(gganimate)

source("get_roc.R")

x <- 1:100
y <- x + rnorm(100, 0, 20)

cutoff <- 80
slider <- 80

point_colour <- function(true, estimated, cutoff, slider) {
    fifelse(true > cutoff & estimated > slider, "True Pos",
            fifelse(true <= cutoff & estimated <= slider, "True Neg",
                    fifelse(true <= cutoff & estimated > slider, "False Pos", "False Neg")))
}




for (slider in seq.int(0L, 100L, by = 10L)) {
    X <- data.table(true = rank(x),
                    estimated = rank(y))

    con_vals <- c("True Pos", "True Neg", "False Pos", "False Neg")

    X[, Top20 := factor(point_colour(true, estimated, cutoff, slider),
                        levels = con_vals)]
    X[, p := factor(fifelse(true > cutoff, 1, 0), levels = c(0, 1))]

    Xroc <- roc(X$p, X$estimated)

    CT <- X[, table(Top20)]
    Se <- CT[["True Pos"]] / (CT[["True Pos"]] + CT[["False Neg"]])
    Sp <- 1 - CT[["False Pos"]] / (CT[["False Pos"]] + CT[["True Neg"]])

    # polys <- data.table(id = as.character(rep(1:4, each = 4)),
    #                     value = rep(1:4, each = 4),
    #                     x = c(0, cutoff, cutoff, 0,
    #                           0, cutoff, cutoff, 0,
    #                           cutoff, 100, 100, cutoff,
    #                           cutoff, 100, 100, cutoff),
    #                     y = c(0, 0, slider, slider,
    #                           slider, slider, 100, 100,
    #                           0, 0, slider, slider,
    #                           slider, slider, 100, 100))

    plt1 <- ggplot(X) +
        # geom_polygon(data = polys, aes(x = x, y = y, fill = id)) +
        geom_point(aes(x = true, y = estimated, colour = Top20), size = 5.0) +
        geom_segment(aes(x = cutoff, xend = cutoff, y = 0, yend = 100)) +
        geom_segment(aes(x = 0, xend = 100, y = slider, yend = slider)) +
        coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) +
        labs(x = "True Rank",
             y = "Estimated Rank",
             title = paste0("Slider = ", slider)) +
        scale_colour_manual(values = c("red", "green", "darkgreen", "darkred"),
                          labels = con_vals,
                          drop = FALSE)
    print(plt1)

    # theme(panel.border = element_rect(colour = "#1b98e0", fill = NA, size = 10))

    plt2 <- ggroc(Xroc) +
        geom_point(aes(x = Xroc$specificities,
                       y = Xroc$sensitivities),
                   size = 0.9,
                   colour = "blue") +
        geom_point(aes(x = Xroc$specificities[slider + 1],
                       y = Xroc$sensitivities[slider + 1]),
                   size = 5.0,
                   colour = "red") +
        labs(x = "1 - Specificity",
             y = "Sensitivity",
             title = paste0("Slider = ", slider))

    plts <- plot_grid(plotlist = list(plt1, plt2), ncol = 2)

    png(filename = paste0("gfx/roc", slider, ".png"), width = 1000, height = 500)
    print(plts)
    dev.off()
}
