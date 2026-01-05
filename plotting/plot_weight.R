library(ggplot2)
library(cowplot)

x <- readRDS("fb_data/fb_12.rds")

wt <- x[sdp == "progeny", .(trial = as.character(trial), weight)]

p1 <- ggplot(wt, aes(x = weight, fill = trial)) +
    geom_density(alpha = 0.5) +
    # scale_x_continuous(trans = "log10") +
    labs(x = "Weight (kg)",
         y = "Density",
         fill = "Trial") +
    theme_classic()
    
p2 <- ggplot(wt, aes(x = log10(weight), fill = trial)) +
    geom_density(alpha = 0.5) +
    labs(x = "log10 Weight (kg)",
         y = "Density",
         fill = "Trial") +
    theme_classic()

plt_wt <- plot_grid(p1, p2, ncol = 1)

ggsave("gfx/weight_density.png",
       plt_wt, width = 6, height = 6)
