library(cowplot)

# These were generated using plot_km_family_trial()
data_t3 <- readRDS("figures/data/fig_km_2.rds")
data_t3 <- readRDS("figures/data/fig_km_8.rds")

scm <- rowwiseDT(
    breaks=, labels=,              values=,
    "fb_d",  "Data Seeder",        "#33A02C",
    "sim_d", "Simulation Seeder",  "#B2DF8A",
    "fb_r",  "Data Contact",       "#1F78B4",
    "sim_r", "Simulation Contact", "#A6CEE3"
)

p1 <- ggplot(data_t3[src == "fb"]) +
    geom_line(aes(x = time, y = survival,
                  group = gp, colour = str,
                  linewidth = src, linetype = src)) +
    scale_fill_manual(breaks = scm$breaks,
                      labels = scm$labels,
                      values = scm$values) +
    scale_colour_manual(breaks = scm$breaks,
                        labels = scm$labels,
                        values = scm$values) +
    scale_linewidth_manual(breaks = c("fb", "sim"),
                           values = c(0.5, 0.2),
                           guide = "none") +
    scale_linetype_manual(breaks = c("fb", "sim"),
                          values = c("solid", "solid"),
                          guide = "none") +
    lims(y = c(0, 1)) +
    # coord_cartesian(expand = FALSE) +
    labs(colour = "Source",
         #linewidth = "Source",
         x = "Time (days)",
         y = "Proportion") +
        facet_grid(cols = vars(variable),
                   rows = vars(trial),
                   scales = "free_x",
                   labeller = labeller(
                       variable = c(Tsign = "(a) Time to visual signs",
                                    RP    = "(b) Time from visual signs to death"),
                       trial = c("1" = "Trial 1",
                                 "2" = "Trial 2"))) +
    theme_bw() +
    theme(legend.position = "bottom",
          strip.background = element_blank())
p1

ggsave("gfx/KM_trials_12.png", p1,
       width = 17.8, height = 17.8 * 9 / 16, units = "cm")
