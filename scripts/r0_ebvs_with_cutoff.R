library(data.table)
library(ggplot2)

x <- fread("fb-final-old/data/scen-2-1-out/ebvs.csv")

setcolorder(x, c(1, 4,5, 2,3, 6,7))
setnames(x, c("id", "sus_g", "sus_e", "inf_g", "inf_e", "tol_g", "tol_e"))
x[, r0 := sus_g + inf_g + tol_g]

q80 <- quantile(x$r0, 0.8) |> unname()

dens <- density(x$r0, adjust = 2)
dd <- with(dens, data.table(x, y))
dd[, Selection := cut(x,
                      breaks = c(min(x), q80, max(x)),
                      labels = c("Lower 80%", "Top 20%"),
                      include.lowest = TRUE)]

gg_colour_hue_r <- function(n) {
    hues <- seq(15, 375, length = n + 1)[- n - 1]
    hcl(h = hues, l = 65, c = 100) |> rev()
}

ggplot(dd, aes(x = x, y = y, fill = Selection)) +
    geom_line() +
    geom_area() +
    geom_vline(aes(xintercept = q80),
               linetype = "dashed") +
    # scale_fill_manual(values = gg_colour_hue_r(2)) +
    labs(x = bquote("EBV for" ~ R[0]),
         y = "Density",
         title = bquote("Density of" ~ R[0] ~ "EBVs, with selection on top 20%")) +
    theme_bw() +
    theme(legend.position = "inside",
          legend.position.inside = c(0.9, 0.8))

ggsave("fb-final/gfx/R0_EBVs.png",
       plot = last_plot(), height = 6, width = 9)

