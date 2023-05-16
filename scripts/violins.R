source("plotting/violins.R")


plts <- plot_violins(data_set = "fb_1_12_Dx_linked", split_replicates = FALSE)
plts_split <- plot_violins(data_set = "fb_1_12_Dx_linked", split_replicates = TRUE)

# plts_split <- plot_violins(data_set = "sim-simple", split_replicates = TRUE)

x <- plts_split
legend <- get_legend(
                     x$beta +
                         guides(color = guide_legend(nrow = 1)) +
                         theme(legend.position = "bottom")
)

plot_set <- c("beta", "latent_period", "detection_period", "recovery_period")

p1 <- plot_grid(plotlist = x[plot_set])
pltg_split <- plot_grid(p1, legend, ncol = 1, rel_heights = c(1, 0.1))
pltg_split

# save(plts_fb_12_Dx_link, plts_fb_12_Dx_split,
#      file="Thursday1.RData")

