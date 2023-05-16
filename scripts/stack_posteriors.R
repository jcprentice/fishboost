source("plotting/posteriors.R")

# Parasites ----
data_set <- "fb-parasites"; num_scenarios <- 12
plts <- list()
for (i in 1:num_scenarios) {
    plts[[i]] <- plot_posteriors(data_set = data_set, scenario = i, ci = "hpdi")
}

# num_scenarios=3 here is Trial 1, Trial 2, Trial 1+2
for (i in 1:3) {
    parms <- names(plts[[i]]$plts)
    splts <- list()
    for (parm in parms) {
        splts[[parm]] <- plot_grid(plts[[i]]$plts[[parm]],
                                   plts[[i+3]]$plts[[parm]],
                                   plts[[i+6]]$plts[[parm]],
                                   plts[[i+9]]$plts[[parm]],
                                   ncol = 1)
        ggsave(glue("gfx/{data_set}/{data_set}-par-{parm}-{i}.pdf"),
               splts[[parm]],
               width = 10, height = 6)
    }
    splts1 <-  plot_grid(plotlist = splts)
    
    ggsave(glue("gfx/{data_set}/{data_set}-stacked_posteriors-{i}.pdf"),
           splts1,
           width = 20, height = 24)
}

ggsave(glue("gfx/{data_set}/{data_set}-par-periods-3.pdf"),
       plot_grid(plotlist = splts[c("beta", "latent_period", "detection_period", "recovery_period")], nrow = 1),
       width = 8, height = 6)

ggsave(glue("gfx/{data_set}/{data_set}-par-trials-3.pdf"),
       plot_grid(plotlist = splts[c("trial_l", "trial_i", "trial_d", "trial_r")], nrow = 1),
       width = 8, height = 6)

ggsave(glue("gfx/{data_set}/{data_set}-par-donors-3.pdf"),
       plot_grid(plotlist = splts[c("donor_l", "donor_i", "donor_d", "donor_r")], nrow = 1),
       width = 8, height = 6)

ggsave(glue("gfx/{data_set}/{data_set}-par-cov_G-3.pdf"),
       plot_grid(plotlist = splts[c("cov_G_ss", "cov_G_ii", "cov_G_rr")], nrow = 1),
       width = 8, height = 6)

ggsave(glue("gfx/{data_set}/{data_set}-par-h2-3.pdf"),
       plot_grid(plotlist = splts[c("h2_ss", "h2_ii", "h2_rr")], nrow = 1),
       width = 8, height = 6)

ggsave(glue("gfx/{data_set}/{data_set}-par-r_G-3.pdf"),
       plot_grid(plotlist = splts[c("r_G_si", "r_G_sr", "r_G_ir")], nrow = 1),
       width = 8, height = 6)

# Censored 1 ----
data_set <- "sim-censored-1"
num_scenarios <- 2
plts <- list()
for (i in 1:num_scenarios) {
    plts[[i]] <- plot_posteriors(data_set = data_set, scenario = i, ci = "hpdi")
}

parms <- names(plts[[1]]$plts)
splts <- list()
for (parm in parms) {
    splts[[parm]] <- plot_grid(plts[[1]]$plts[[parm]],
                               plts[[2]]$plts[[parm]],
                               # plts[[3]]$plts[[parm]],
                               ncol = 1)
}
splts1 <- plot_grid(plotlist = splts)

ggsave(glue("gfx/{data_set}/{data_set}-stacked_posteriors.pdf"),
       splts1,
       width = 10, height = 12)


