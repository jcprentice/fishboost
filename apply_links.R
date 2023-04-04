# Take params and ensure that linked traits are correct by setting Sigma_G
apply_links <- function(params) {
    {
        all_traitnames  <- params$all_traitnames
        traitnames      <- params$traitnames
        ntraits         <- params$ntraits
        sim_trial_fe    <- params$sim_trial_fe
        sim_donor_fe    <- params$sim_donor_fe
        sim_txd_fe      <- params$sim_txd_fe
        trial_fe        <- params$trial_fe
        donor_fe        <- params$donor_fe
        txd_fe          <- params$txd_fe
        sim_link_traits <- params$sim_link_traits
        sim_link_trial  <- params$sim_link_trial
        sim_link_donor  <- params$sim_link_donor
        sim_link_txd    <- params$sim_link_txd
        sim_link_shapes <- params$sim_link_shapes
        link_traits     <- params$link_traits
        link_trial      <- params$link_trial
        link_donor      <- params$link_donor
        link_txd        <- params$link_txd
        link_shapes     <- params$link_shapes
        fe_vals         <- params$fe_vals
        Sigma_G         <- params$Sigma_G
        Sigma_E         <- params$Sigma_E
        r_eta_shape     <- params$r_eta_shape
        r_rho_shape     <- params$r_rho_shape
        r_gamma_shape   <- params$r_gamma_shape
        eta_type        <- params$eta_type
        rho_type        <- params$rho_type
        gamma_type      <- params$gamma_type
    }
    
    # Need to make explicit copies because datatables work by reference
    params2 <- copy(params)
    
    # All changes made to priors also happen in params2
    priors  <- params2$priors
    
    # Grab the first letter of each traitname for convenience
    traitnames1 <- substr(traitnames, 1, 1)
    all_traits1 <- substr(all_traitnames, 1, 1)
    
    
    # Individual Effects ----
    if (ntraits > 0) {
        # Turn sim_link_traits into string vector
        z <- strsplit(sim_link_traits, "")[[1]]
        
        # Subset based on which traits are actually used
        # z <- z[match(traitnames1, all_traits1)]
        
        # Turn that into indexes for Sigmas
        # new_idxs <- match(z, traitnames1)
        new_idxs <- match(z, all_traits1)
        
        # Copy rows and columns in Sigmas
        Sigma_G <- Sigma_G[new_idxs, , drop = FALSE]
        Sigma_G <- Sigma_G[, new_idxs, drop = FALSE]
        Sigma_E <- Sigma_E[new_idxs, , drop = FALSE]
        Sigma_E <- Sigma_E[, new_idxs, drop = FALSE]
        
        # Fix values in priors
        cor_G <- cov2cor(Sigma_G)
        cor_E <- cov2cor(Sigma_E)
        
        iS <- which(traitnames1 == "s")
        iI <- which(traitnames1 == "i")
        iR <- which(traitnames1 == "r")
        
        # Be careful here in case any traits aren't used
        priors[parameter == "cov_G_ss", true_val := if (use) Sigma_G[iS, iS] else true_val]
        priors[parameter == "cov_G_ii", true_val := if (use) Sigma_G[iI, iI] else true_val]
        priors[parameter == "cov_G_rr", true_val := if (use) Sigma_G[iR, iR] else true_val]
        priors[parameter == "r_G_si",   true_val := if (use) cor_G[iS, iI] else true_val]
        priors[parameter == "r_G_sr",   true_val := if (use) cor_G[iS, iR] else true_val]
        priors[parameter == "r_G_ir",   true_val := if (use) cor_G[iI, iR] else true_val]
        priors[parameter == "cov_E_ss", true_val := if (use) Sigma_E[iS, iS] else true_val]
        priors[parameter == "cov_E_ii", true_val := if (use) Sigma_E[iI, iI] else true_val]
        priors[parameter == "cov_E_rr", true_val := if (use) Sigma_E[iR, iR] else true_val]
        priors[parameter == "r_E_si",   true_val := if (use) cor_E[iS, iI] else true_val]
        priors[parameter == "r_E_sr",   true_val := if (use) cor_E[iS, iR] else true_val]
        priors[parameter == "r_E_ir",   true_val := if (use) cor_E[iI, iR] else true_val]
        
        # Copy updated values back into params2
        params2$Sigma_G <- Sigma_G
        params2$Sigma_E <- Sigma_E
        params2$h2 <- Sigma_G / (Sigma_G + Sigma_E)
        params2$cor_G <- cor_G
        params2$cor_E <- cor_E
    }
    
    
    # Fixed Effects ----
    if (sim_trial_fe != "") {
        # Copy FE values across, overwriting previous value
        sim_link_trial_v = strsplit(sim_link_trial, "")[[1]]
        zt <- match(sim_link_trial_v, all_traits1)
        fe_vals["trial", ] <- fe_vals["trial", zt]
    }
    
    if (trial_fe != "") {
        # sieve SLIDR for trial_fe, perform the map, take unique values and
        # convert into list of trial priors
        # e.g. "ld" -> ".l.d." -> ".l.l" -> "ll" -> "l" -> "trial_l"
        link_trial_v = strsplit(link_trial, "")[[1]]
        trial_priors <- paste0("trial_", unique(link_trial_v[grepl(paste0("[", trial_fe, "]"), all_traits1)]))
        priors[startsWith(parameter, "trial_"), use := parameter %in% trial_priors]
    }
    
    if (sim_donor_fe != "") {
        sim_link_donor_v = strsplit(sim_link_donor, "")[[1]]
        zd <- match(sim_link_donor_v, all_traits1)
        fe_vals["donor", ] <- fe_vals["donor", zd]
    }
    
    if (donor_fe != "") {
        link_donor_v = strsplit(link_donor, "")[[1]]
        donor_priors <- paste0("donor_", unique(link_donor_v[grepl(paste0("[", donor_fe, "]"), all_traits1)]))
        priors[startsWith(parameter, "donor_"), use := parameter %in% donor_priors]
    }
    
    if (sim_txd_fe != "") {
        sim_link_txd_v = strsplit(sim_link_txd, "")[[1]]
        zd <- match(sim_link_txd_v, all_traits1)
        fe_vals["txd", ] <- fe_vals["txd", zd]
    }
    
    if (txd_fe != "") {
        link_txd_v = strsplit(link_txd, "")[[1]]
        txd_priors <- paste0("txd_", unique(link_txd_v[grepl(paste0("[", txd_fe, "]"), all_traits1)]))
        priors[startsWith(parameter, "txd_"), use := parameter %in% txd_priors]
    }
    
    # Copy updated values back into params2
    params2$fe_vals <- fe_vals
    
    # Shape parameters ----
    shapes <- c(l = r_eta_shape, d = r_rho_shape, r = r_gamma_shape)
    shape <- c(lat = shapes[[substr(sim_link_shapes, 1, 1)]],
               det = shapes[[substr(sim_link_shapes, 2, 2)]],
               rec = shapes[[substr(sim_link_shapes, 3, 3)]])
    
    # Copy updated values back into params2
    params2$r_eta_shape <- shape[["lat"]]
    params2$r_rho_shape <- shape[["det"]]
    params2$r_gamma_shape <- shape[["rec"]]
    
    # Set values in priors (be careful not to accidentally turn them back on if the type is "exp")
    priors[parameter == "eta_shape",   use := use && grepl("l", link_shapes)]
    priors[parameter == "rho_shape",   use := use && grepl("d", link_shapes)]
    priors[parameter == "gamma_shape", use := use && grepl("r", link_shapes)]
    
    return(params2)
}
