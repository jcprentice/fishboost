rename_pars <- function(pars) {
    pars <- sub("latent_period", "Latent Period (days)", pars)
    pars <- sub("eta_shape", "LP shape", pars)
    pars <- sub("detection_period", "Detection Period (days)", pars)
    pars <- sub("rho_shape", "DP shape", pars)
    pars <- sub("recovery_period", "Recovery Period (days)", pars)
    pars <- sub("gamma_shape", "RP shape", pars)
    pars <- sub("_l$", "_Latency", pars)
    pars <- sub("_i$", "_Infectivity", pars)
    pars <- sub("_d$", "_Detectability", pars)
    pars <- sub("_r$", "_Recoverability", pars)
    pars <- sub("trial_", "Trial ", pars)
    pars <- sub("donor_", "Donor ", pars)
    pars <- sub("sigma", "Group Effect", pars)
    pars <- sub("^cov_", "Var_", pars)
    pars <- sub("^r_", "Cor_", pars)
    pars <- sub("_ss$", "(Sus)", pars)
    pars <- sub("_ii$", "(Inf)", pars)
    pars <- sub("_rr$", "(Rec)", pars)
    pars <- sub("_si$", "(Sus, Inf)", pars)
    pars <- sub("_sr$", "(Sus, Rec)", pars)
    pars <- sub("_ir$", "(Inf, Rec)", pars)
    # params2 <- gsub("_", " ", params2)
    
    pars
}

param_order <- c("beta",
                 "latent_period", "eta_shape",
                 "detection_period", "rho_shape",
                 "recovery_period", "gamma_shape",
                 "sigma",
                 "cov_G_ss", "cov_G_ii", "cov_G_rr", "r_G_si", "r_G_sr", "r_G_ir",
                 "cov_E_ss", "cov_E_ii", "cov_E_rr", "r_E_si", "r_E_sr", "r_E_ir",
                 "trial_l", "trial_i",  "trial_d", "trial_r",
                 "donor_l", "donor_i",  "donor_d", "donor_r")