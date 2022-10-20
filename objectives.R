# library(ggplot2)
library(tidyverse)
library(ggthemes)

{
    obj <- read.csv("objectives.csv")

    # delete empty rows
    obj <- obj[!is.na(obj$Objective), ]

    # tidy up and format
    obj$Task      <- factor(obj$Task, levels = rev(obj$Task)) # ordered = TRUE ?
    obj$Start     <- as.Date(obj$Start)
    obj$End       <- as.Date(obj$End)
    obj$Objective <- factor(obj$Objective)
}

gantt_out <- ggplot(obj,
                    aes(x = Start, xend = End,
                        y = Task, yend = Task,
                        colour = Objective)) +
    geom_segment(size = 8) +
    geom_vline(xintercept = Sys.Date(),
               linetype = "dashed") +
    labs(title = "Jamie Prentice PDR Objectives",
         x = "Time", y = "Task") +
    theme(#axis.text = element_text(size = 12),
          legend.position = "bottom")

gantt_out
