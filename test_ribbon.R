library(tidyverse)

tbl <- tibble(
    x = seq(-2, 2, by = 1e-2),
    y = dnorm(x, 0, 0.5),
    q = factor(cut(x, qnorm(seq(0, 1, by = 0.1), 0, 0.5),
                   right = FALSE, labels = FALSE, include.lowest = TRUE)))

ggplot(tbl) +
    geom_ribbon(aes(x = x, ymin = 0, ymax = y, fill = q)) +
    scale_fill_viridis_d()
