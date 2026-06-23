library(ggplot2)
library(ggtext)

theme_jamie <- function() {
    theme_classic() +
    theme(legend.position = "none",
          text = element_text(size = 5),
          plot.title = element_markdown(size = 7),
          # panel.background = element_blank(),
          # panel.grid.major = element_blank(),
          # panel.grid.minor = element_blank(),
          axis.line = element_line(linewidth = 1 / .pt),
          axis.ticks = element_line(linewidth = 0.5 / .pt))
}
