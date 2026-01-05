library(data.table)
library(ggplot2)
library(cowplot)

x <- readRDS("fb_data/fb_12.rds")
x[, trial := factor(trial)]

p1 <- ggplot(x) +
    geom_density(aes(x = weight, fill = trial),
                   position = "identity", alpha = 0.5) +
    theme(legend.position = c(0.9, 0.8))
p2 <- ggplot(x) +
    geom_density(aes(x = log(weight), fill = trial),
                   position = "identity", alpha = 0.5) +
    theme(legend.position = c(0.9, 0.8))


p1 <- ggplot(x) +
    geom_histogram(aes(x = weight, y = after_stat(density), fill = trial),
                   position = "identity", alpha = 0.5, bins = 100) +
    theme(legend.position = c(0.9, 0.9))

p2 <- ggplot(x) +
    geom_histogram(aes(x = log10(weight), y = after_stat(density), fill = trial),
                   position = "identity", alpha = 0.5, bins = 100) +
    theme(legend.position = c(0.9, 0.9))

plt <- plot_grid(p1, p2, ncol = 1)
ggsave("weight_density.png", plt, width = 6, height = 5)

d1 = rnorm(1000, mean = 50, sd = 10)
d2 = 2 * d1
x = data.table(trial = rep(factor(c(1,2)), each = 1000),
               weight = c(d1, d2))
