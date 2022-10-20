library(fitdistrplus)

fb_traits <- readRDS("fb_data/fb_traits.rds")

# Latent Periods
LP  <- fb_traits[!is.na(Tinf) & !is.na(Tsym),              Tsym - Tinf]
LP1 <- fb_traits[!is.na(Tinf) & !is.na(Tsym) & trial == 1, Tsym - Tinf]
LP2 <- fb_traits[!is.na(Tinf) & !is.na(Tsym) & trial == 2, Tsym - Tinf]

# Recovery Periods
RP1 <- fb_traits[!is.na(Tsym) & !is.na(Trec) & Trec < 104 & trial == 1, Trec - Tsym]
RP2 <- fb_traits[!is.na(Tsym) & !is.na(Trec) & Trec < 160 & trial == 2, Trec - Tsym]
RP <- c(RP1, RP2)

# Overall
fitdist(LP + 0.5, "gamma")
fitdist(RP + 0.5, "gamma")

# Just Trial 1
LP1fit <- fitdist(LP1 + 0.5, "gamma")
RP1fit <- fitdist(RP1 + 0.5, "gamma")

# These are the values we'll actually use
r_eta_shape <- LP1fit$estimate[["shape"]]
r_eta_rate  <- LP1fit$estimate[["rate"]]
r_eta  <- r_eta_rate / r_eta_shape
mean_LP <- 1 / r_eta

print(signif(c(LP = mean_LP, shape = r_eta_shape, rate = r_eta_rate), 3))

r_gamma_shape <- RP1fit$estimate[["shape"]]
r_gamma_rate  <- RP1fit$estimate[["rate"]]
r_gamma  <- r_gamma_rate / r_gamma_shape
mean_RP <- 1 / r_gamma

print(signif(c(RP = mean_RP, shape = r_gamma_shape, rate = r_gamma_rate), 3))

# Just Trial 2
fitdist(LP2 + 0.5, "gamma")
fitdist(RP2 + 0.5, "gamma")
