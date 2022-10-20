library(data.table)
library(ggplot2)
library(cowplot)

source("make_parameters.R")
source("make_time_series.R")
source("make_plots.R")
source("add_latency.R")

fb_traits <- readRDS("fb_data/fb_traits.rds")

# Work with copy of data
fb <- fb_traits[sdp == "progeny", .(id, sdp, trial, group, donor, Tinf, Tsym, Trec)]
fb[Tsym > 10, donor := 0]
fb[, Tinf := fifelse(donor == 1, 0, NA_real_)]


# Patch in missing values
# add_latency_donors(fb)
add_latency_donors_plus(fb)
fb[, `:=`(Tinf = estimated_Tinf, estimated_Tinf = NULL)]
fb[group == 3]

# Plot as though SEIR model
{
    params <- make_parameters(model = "SIDR", use_fb_data = TRUE)
    plt_fb1 <- plot_model(fb[trial == 1], params)
    plt_fb2 <- plot_model(fb[trial == 2], params)
    plt_fb <- plot_model(fb, params)
    print(plt_fb1)
}

plot_paths <- function(pop, group_num = 3, sort_ids = TRUE) {
    pop1 <- pop[group %in% group_num]

    if ("estimated_Tinf" %in% names(pop1)) {
        setnames(pop1, "estimated_Tinf", "Tinf")
    }

    if (sort_ids) {
        pop1[order(-donor, Tsym, Trec), id := sort(id)]
    }

    plt1 <- ggplot(pop1) +
        geom_segment(aes(x = Tsym, xend = Trec, y = id, yend = id, colour = "Symptomatic"), size = 1) +
        geom_point(aes(x = Tsym, y = id, colour = "Tsym"), size = 3) +
        geom_point(aes(x = Trec, y = id, colour = "Trec"), size = 3) +
        scale_y_reverse() +
        labs(x = "Time (days)",
             title = paste0("FB group ", group_num, " with symptom and removal times")) +
        scale_colour_manual("Compartments",
                            breaks = c("Symptomatic", "Tsym", "Trec"),
                            values = c("blue", "blue", "purple"),
                            guide = guide_legend(override.aes = list(
                                linetype = c(rep("solid", 1), rep("blank", 2)),
                                shape = c(rep(NA, 1), rep(16, 2))))) #+
    # theme(legend.position = "bottom")

    plt2 <- ggplot(pop1) +
        geom_segment(aes(x = Tinf, xend = Trec, y = id, yend = id, colour = "Overall"), linetype = 2) +
        geom_segment(aes(x = Tinf, xend = Tsym, y = id, yend = id, colour = "Hidden"), size = 1) +
        geom_segment(aes(x = Tsym, xend = Trec, y = id, yend = id, colour = "Symptomatic"), size = 1) +
        geom_point(aes(x = Tinf, y = id, colour = "Estimated Tinf"), size = 3) +
        geom_point(aes(x = Tsym, y = id, colour = "Tsym"), size = 3) +
        geom_point(aes(x = Trec, y = id, colour = "Trec"), size = 3) +
        scale_y_reverse() +
        labs(x = "Time (days)",
             title = paste0("FB group ", group_num, " with inferred infection times")) +
        scale_colour_manual("Compartments",
                            breaks = c("Overall", "Hidden", "Symptomatic", "Estimated Tinf", "Tsym", "Trec"),
                            values = c("black", "red", "blue", "red", "blue", "purple"),
                            guide = guide_legend(override.aes = list(
                                linetype = c("dashed", rep("solid", 2), rep("blank", 3)),
                                shape = c(rep(NA, 3), rep(16, 3))))) #+
    # theme(legend.position = "bottom")

    plt <- plot_grid(plotlist = list(plt1, plt2))

    print(plt)

    list(plt1, plt2, plt)
}

{
    x <- plot_paths(fb, group_num = 3, sort_ids = FALSE)

    png("gfx/group3-noLP.png", width = 1000, height = 1000)
    print(x[1])
    dev.off()

    png("gfx/group3-LP.png", width = 1000, height = 1000)
    print(x[2])
    dev.off()

    png("gfx/group3.png", width = 2000, height = 1000)
    print(x[3])
    dev.off()

    x <- plot_paths(fb, group_num = 3, sort_ids = TRUE)

    png("gfx/group3-noLP-sorted.png", width = 1000, height = 1000)
    print(x[1])
    dev.off()

    png("gfx/group3-LP-sorted.png", width = 1000, height = 1000)
    print(x[2])
    dev.off()

    png("gfx/group3-sorted.png", width = 2000, height = 1000)
    print(x[3])
    dev.off()
}


 plot_grid(plotlist=list(plt_fb1, plt))
