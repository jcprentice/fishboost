library(glue)

patch_params <- function(params, trace_row = 0) {
    data_set        <- params$patch_data_set
    scenario        <- params$patch_scenario
    patch_with_mean <- params$patch_with_mean
    
    if (is.na(data_set) || data_set == "") {
        message("No patches necessary")
        return(params)
    }
    
    # Make changes to a copy of params and return that
    params2 <- copy(params)
    
    message(glue("Patching params with {ptype} posteriors from '{data_set}' scenario {scenario}...",
                 ptype = if (patch_with_mean) "mean" else "sampled"))
    
    f <- glue("data/{data_set}/scen-{scenario}-1_out/trace_combine.txt")
    if (!file.exists(f)) {
        message("- missing trace file")
        return(params)
    }
    trace_file <- fread(f)
    trace_file[, names(trace_file)[grepl("Group effect", names(trace_file))] := NULL]
    
    # Extract list of posterior values from trace file
    posteriors <- if (patch_with_mean) {
        lapply(trace_file, mean)[-1]
    } else {
        if (trace_row <= 0) {
            trace_row = sample(nrow(trace_file), 1)
        }
        message(glue(" - using row {trace_row} from trace file"))
        as.list(trace_file[trace_row])[-1]
    }
    rm(trace_file)
    
    
    # Load results file to get any parameters not saved in trace file
    f <- glue("results/{data_set}/scen-{scenario}-1.RData")
    results_exists <- file.exists(f)
    if (results_exists) {
        message(" - extracting parameters from results file")
        tmp_env <- new.env()
        load(f, envir = tmp_env)
        priors <- tmp_env$params$priors
    } else {
        message(" - no results file")
        priors <- params$priors
    }
    
    # Copy prior values, using true_val from the posteriors
    params2$priors[, `:=`(type = priors$type,
                          val1 = priors$val1,
                          val2 = priors$val2,
                          true_val = unlist(posteriors)[parameter],
                          use = priors$use)]
    
    # Patch from Trace file
    
    pars_patched <- "priors"
    
    sir <- c(s = "susceptibility", i = "infectivity", r = "recoverability")
    patched_covariances <- FALSE
    for (GE in c("G", "E")) for (x in names(sir)) for (y in names(sir)) {
        r_xy <- glue("{rcov}_{GE}_{x}{y}",
                     rcov = ifelse(x == y, "cov", "r"))
        cov_xx <- glue("cov_{GE}_{x}{x}")
        cov_yy <- glue("cov_{GE}_{y}{y}")

        if (!r_xy %in% names(posteriors))  next

        cov_xy <- if (x == y) {
            with(posteriors, get(cov_xx))
        } else {
            with(posteriors, get(r_xy) * sqrt(get(cov_xx) * get(cov_yy)))
        }

        params2[[glue("Sigma_{GE}")]][sir[[x]], sir[[y]]] <- cov_xy
        params2[[glue("Sigma_{GE}")]][sir[[y]], sir[[x]]] <- cov_xy
        patched_covariances <- TRUE
    }
    if (patched_covariances) {
        pars_patched <- c(pars_patched, "covariances")
    }
    
    
    # if (all(c("r_G_si", "r_G_sr", "r_G_ir") %in% names(posteriors))) {
    #     pars_patched <- c(pars_patched, "covariances")
    #     
    #     # for indexing Sigma
    #     sir <- c("susceptibility", "infectivity", "recoverability")
    #     
    #     cors_G <- matrix(with(posteriors, c(1, r_G_si, r_G_sr, r_G_si, 1, r_G_ir, r_G_sr, r_G_ir, 1)), 3, 3)
    #     covs_G <- diag(with(posteriors, sqrt(c(cov_G_ss, cov_G_ii, cov_G_rr))), 3, 3)
    #     Sigma_G <- covs_G %*% cors_G %*% covs_G
    #     dimnames(Sigma_G) <- list(sir, sir)
    #     
    #     cors_E <- matrix(with(posteriors, c(1, r_E_si, r_E_sr, r_E_si, 1, r_E_ir, r_E_sr, r_E_ir, 1)), 3, 3)
    #     covs_E <- diag(with(posteriors, sqrt(c(cov_E_ss, cov_E_ii, cov_E_rr))), 3, 3)
    #     Sigma_E <- covs_E %*% cors_E %*% covs_E
    #     dimnames(Sigma_E) <- list(sir, sir)
    #     
    #     params$Sigma_G[sir, sir] <- Sigma_G
    #     params$Sigma_E[sir, sir] <- Sigma_E
    #     
    #     for (p in c("r_G_si", "r_G_sr", "r_G_si", "r_G_ir", "r_G_sr", "r_G_ir",
    #                 "r_G_si", "r_G_sr", "r_G_si", "r_G_ir", "r_G_sr", "r_G_ir",
    #                 "cov_G_ss", "cov_G_ii", "cov_G_rr",
    #                 "cov_E_ss", "cov_E_ii", "cov_E_rr")) {
    #         params$priors[parameter == "p", true_val := get(p)]
    #     }
    # }
    
    
    
    
    if ("beta" %in% names(posteriors)) {
        params2$r_beta <- posteriors$beta
        pars_patched <- c(pars_patched, "beta")
    }
    if ("latent_period" %in% names(posteriors)) {
        params2$r_eta <- params2$r_eta_rate <- 1 / posteriors$latent_period
        pars_patched <- c(pars_patched, "latent period")
    }
    if ("eta_shape" %in% names(posteriors)) {
        params2$r_eta_shape <- posteriors$eta_shape
        params2$r_eta_rate <- params2$r_eta_rate * posteriors$eta_shape
        pars_patched <- c(pars_patched, "eta shape")
    }
    if ("detection_period" %in% names(posteriors)) {
        params2$r_rho <- params2$r_rho_rate <- 1 / posteriors$detection_period
        pars_patched <- c(pars_patched, "detection period")
    }
    if ("rho_shape" %in% names(posteriors)) {
        params2$r_rho_shape <- posteriors$rho_shape
        params2$r_rho_rate <- params2$r_rho_rate * posteriors$rho_shape
        pars_patched <- c(pars_patched, "rho shape")
    }
    if ("recovery_period" %in% names(posteriors)) {
        params2$r_gamma <- params2$r_gamma_rate <- 1 / posteriors$recovery_period
        pars_patched <- c(pars_patched, "recovery period")
    }
    if ("gamma_shape" %in% names(posteriors)) {
        params2$r_gamma_shape <- posteriors$gamma_shape
        params2$r_gamma_rate <- params2$r_gamma_rate * posteriors$gamma_shape
        pars_patched <- c(pars_patched, "gamma shape")
    }
    if ("sigma" %in% names(posteriors)) {
        # params2$group_effect <- posteriors$sigma
        pars_patched <- c(pars_patched, "group effect")
    }
    
    
    for (fe_type in c("trial", "donor", "txd")) {
        for (fe_trait in params$all_traitnames) {
            # build name e.g. trial_l
            fe_name <- paste0(fe_type, "_", substr(fe_trait, 1, 1))
            if (fe_name %in% names(posteriors)) {
                params2$fe_vals[fe_type, fe_trait] <- posteriors[[fe_name]]
                pars_patched <- c(pars_patched, fe_name)
            }
        }
    }
    
    if (results_exists) {
        # Ensure FEs match
        params2$trial_fe <- tmp_env$params$trial_fe
        params2$donor_fe <- tmp_env$params$donor_fe
        params2$txd_fe   <- tmp_env$params$txd_fe
        pars_patched    <- c(pars_patched, "fixed effects")
        
        # shape parameters
        params2$r_eta_shape <- tmp_env$params$r_eta_shape
        params2$r_rho_shape <- tmp_env$params$r_rho_shape
        params2$r_gamma_shape <- tmp_env$params$r_gamma_shape
        
        # I think this is necessary?
        params2$time_step <- tmp_env$params$time_step
        params2$censor <- tmp_env$params$censor
        pars_patched <- c(pars_patched, "time_step", "censor")
    }
    
    # Check to make sure priors are in bounds, don't bother with sigma
    tmp <- as.list(params2$priors[parameter == "sigma"])
    params2$priors[val1 >= true_val, val1 := true_val - 2]
    params2$priors[val2 <= true_val, val2 := true_val + 2]
    params2$priors[parameter == "sigma", (colnames(params$priors)) := tmp]
    
    message(" - ", paste(pars_patched, collapse = ", "))
    
    params2
}
