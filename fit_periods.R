library(fitdistrplus)

fb_data <- readRDS("fb_data/fb_data12.rds")

{
    x <- fb_data[sdp == "progeny", .(trial, Tsym, Trec)]
    
    x[, Trec2 := Trec]
    x[Trec2 %in% c(104, 160), Trec2 := NA]
    
    # Cut out values where both Tsym and Trec2 are NA
    x <- x[!is.na(Tsym) | !is.na(Trec2)]
    
    minRP <- 0.1
    
    # Put in naive values
    x[, `:=`(left = pmax(Trec - Tsym, minRP),
             right = Trec2 - Tsym)]
    
    # Fix right
    # Note: Trec2 = NA means Trec2 > 104 or 160
    x[is.na(Tsym) & !is.na(Trec2), right := Trec2]
    # x[is.na(Tsym) & !is.na(Trec2), right := minRP]
    
    # Fix left
    x[is.na(Tsym) & !is.na(Trec2), left := Trec2 - minRP]
    # x[is.na(Tsym) & !is.na(Trec2), left := minRP]
    
    # Fix zero RP
    x[Tsym == Trec, `:=`(left = minRP, right = minRP)]
    
    
    print(x[1:25, .(Tsym, Trec2, left, right, delta = right-left)])
    
    
    f1 <- fitdistcens(x[trial == 1] + 0.5, "gamma")
    f1$estimate[["mean"]] = f1$estimate[["shape"]] / f1$estimate[["rate"]]
    message("Trial 1:")
    print(round(f1$estimate, 1))
    
    f2 <- fitdistcens(x[trial == 2] + 0.5, "gamma")
    f2$estimate[["mean"]] = f2$estimate[["shape"]] / f2$estimate[["rate"]]
    message("Trial 2:")
    print(round(f2$estimate, 1))
    
    f12 <- fitdistcens(x + 0.5, "gamma")
    f12$estimate[["mean"]] = f12$estimate[["shape"]] / f12$estimate[["rate"]]
    message("Trial 1+2:")
    print(round(f12$estimate, 1))
}



# Latent Periods
LP  <- fb_data[!is.na(Tinf) & !is.na(Tsym),              Tsym - Tinf]
LP1 <- fb_data[!is.na(Tinf) & !is.na(Tsym) & trial == 1, Tsym - Tinf]
LP2 <- fb_data[!is.na(Tinf) & !is.na(Tsym) & trial == 2, Tsym - Tinf]

# Recovery Periods
RP1 <- fb_data[!is.na(Tsym) & !is.na(Trec) & Trec < 104 & trial == 1, Trec - Tsym]
RP2 <- fb_data[!is.na(Tsym) & !is.na(Trec) & Trec < 160 & trial == 2, Trec - Tsym]
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
