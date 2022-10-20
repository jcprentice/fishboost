library(data.table)
library(ggplot2)
library(cowplot)

source("make_parameters.R")
source("make_time_series.R")
source("make_plots.R")
source("add_latency.R")

fb_traits <- readRDS("fb_data/fb_traits.rds")

# Work with copy of data ----
fb <- fb_traits[sdp == "progeny", .(id, sdp, trial, group, donor, Tinf, Tsym, Trec)]
fb[, Tinf := fifelse(donor == 1, 0, NA_real_)]


# Patch in missing values ----

# add_latency_donors(fb)
# add_latency_donors_plus(fb)
if ("estimated_Tinf" %in% names(fb)) {
    fb[, `:=`(Tinf = estimated_Tinf, estimated_Tinf = NULL)]
}
# fb[group == 3]

dt_donors <- fb[!is.na(Tsym) & donor == 1, .(Trial = as.factor(trial), Tsym)]
dt_all <- fb[!is.na(Tsym), .(Trial = as.factor(trial), Tsym)]

# manual Y scaling
d_donors <- fb[donor == 1 & !is.na(Tsym), .N, trial][, N] / 180^2
d_all    <- fb[!is.na(Tsym), .N, trial][, N] / 900^2

my_breaks <- seq.int(0L, 154L, by = 1L)
n_breaks  <- length(my_breaks) - 1L

scale_donors <- rep(d_donors, each = n_breaks)
scale_all    <- rep(d_all,    each = n_breaks)


plt_donors <- ggplot(data = dt_donors, aes(x = Tsym, fill = Trial, colour = Trial)) +
    geom_histogram(aes(y = scale_donors * ..count..),
                   alpha = 0.5,
                   position = "identity",
                   # breaks = c(0:10, seq.int(20L, 160L, by = 10L)))
                   breaks = my_breaks) +
    # geom_density(alpha = 0., adjust = 1) +
    labs(title = "Donor Tsyms for Trials 1 & 2",
         x = "Time until symptoms (days)",
         y = "Density")

plt_all <- ggplot(data = dt_all, aes(x = Tsym, fill = Trial, colour = Trial)) +
    geom_histogram(aes(y = scale_all * ..count..),
                   alpha = 0.5,
                   position = "identity",
                   # breaks = c(0:10, seq.int(20L, 160L, by = 10L)))
                   breaks = my_breaks) +
    # geom_density(alpha = 0., adjust = 1) +
    labs(title = "All Tsyms for Trials 1 & 2",
         x = "Time until symptoms (days)",
         y = "Density")

plot_grid(plotlist = list(plt_donors, plt_all), ncol=1)
