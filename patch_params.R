{
    library(data.table)
    library(stringr)
    library(purrr)
    
    source("widen_priors.R")
    source("remove_bici_fes.R")
}

# This will read an existing fitted dataset and apply either the mean or median
# parameter values, or sample from the posteriors. It saves the source so that
# the corresponding BVs can be applied too (but we do that later).

patch_params <- function(params, trace_row = 0) {
    {
        patch_dataset <- params$patch_dataset %||% ""
        patch_name    <- params$patch_name
        patch_type    <- params$patch_type %||% "median"
        patch_state   <- params$patch_state %||% FALSE
        model_traits  <- params$model_traits
        use_traits    <- params$use_traits
        skip_patches  <- params$skip_patches
        sim_new_data  <- params$sim_new_data
        msgs          <- params$msgs
    }

    if (msgs) message("Patching parameters ...")

    if (is.na(patch_dataset) || patch_dataset == "") {
        if (msgs) message("- no patches necessary")
        return(params)
    }
    
    if (msgs) message(str_glue("- using *{patch_type}* posteriors from ",
                               "'{patch_dataset}' scenario {patch_name} ..."))
    
    {
        data_dir    <- str_glue("datasets/{patch_dataset}/data")
        results_dir <- str_glue("datasets/{patch_dataset}/results")
        out_dir     <- str_glue("{data_dir}/{patch_name}-out")
        inf_out_dir <- str_glue("{out_dir}/output-inf")
        sim_out_dir <- str_glue("{out_dir}/output-sim")
    }
    
    # Make changes to a copy of params and later return that
    params2 <- copy(params)
    
    
    # If 'patch_state' is specified, use that, otherwise check if a trace file
    # is present and store the row used.
    
    etc_src <- if (sim_new_data == "etc_sim") "etc_sim" else "etc_inf"
    
    # Select patch ----
    
    # This might be already set, in which case skip
    if (!is.numeric(patch_state) && patch_state == TRUE) {
        f <- str_glue("{out_dir}/{etc_src}.rds")
        patch_state <- if (!file.exists(f)) {
            message("no state files found, continuing without")
            FALSE
        } else {
            if (patch_type == "sampled") {
                readRDS(f)$parameters[, sample(state, 1)]
            } else {
                0
            }
        }
    }
    params2$patch_state <- patch_state
    
    # Load patch ----
    
    # Use patch_state
    if (is.numeric(patch_state)) {
        if (msgs) message(str_glue("- using state: {patch_state}"))
        
        f <- str_glue("{out_dir}/{etc_src}.rds")
        tmp <- readRDS(f)$parameters[!str_starts(parameter, "^Group")]
        tmp <- if (patch_state == 0) {
            tmp[, .(value = get(patch_type)(value)), parameter]
        } else {
            tmp[state == patch_state]
        }
        # This is to under renaming the parameters in BICI
        patch_vals <- setNames(as.list(tmp$value),
                               tmp$parameter |> str_replace_all(
                                   c("LP" = "latent_period",
                                     "DP" = "detection_period",
                                     "RP" = "removal_period")))
        rm(tmp)
    }
    
    # Use trace file
    # Caution: 0 == FALSE, so instead test if is logical, since we've already
    # handled is TRUE
    if (is.logical(patch_state)) {
        # Get parameters from trace file
        f <- str_glue("{out_dir}/extended_trace_combine.tsv")
        if (!file.exists(f)) {
            f <- str_remove(f, "extended_")
            if (!file.exists(f)) {
                stop("No trace file")
            }
        }
        if (msgs) message(str_glue("- reading trace file '{f}'"))
        trace_file <- fread(f)
        # trace_file[, str_subset(names(trace_file), "^\\d", negate = TRUE) := NULL]
        
        # Extract list of posterior values from trace file
        patch_vals <- if (patch_type == "mean") {
            map(trace_file, mean)
        } else if (patch_type == "median") {
            map(trace_file, median)
        } else {
            if (!trace_row %between% c(1, nrow(trace_file))) {
                trace_row <- sample(nrow(trace_file), 1)
            }
            params2$trace_row <- trace_row
            if (msgs) message(str_glue("- using row {trace_row} from trace file"))
            as.list(trace_file[trace_row])
        }
        rm(trace_file)
    }
    
    
    # Load results file to get any parameters not saved in trace file
    rf <- str_glue("{results_dir}/{patch_name}.rds")
    # Save this for later
    resfile_exists <- file.exists(rf)
    if (resfile_exists) {
        if (msgs) message(str_glue("- extracting parameters from '{rf}'"))
        tmp_params <- readRDS(rf)$params
        priors <- tmp_params$priors
    } else {
        if (msgs) message("- no results file")
        priors <- params$priors
    }
    
    # Correct for old type
    priors[type == "Flat",  type := "uniform"]
    priors[type == "Fixed", type := "constant"]
    
    # Determine what to patch ----
    
    # Copy prior values, using `true_val` from the patch_vals
    
    pp <- intersect(setdiff(names(patch_vals), skip_patches),
                    params2$priors$parameter)
    
    # Now remove anything we don't want to patch
    walk(skip_patches, \(x) {
        pp <<- pp |>
            str_subset(switch(x,
                              "base" = "sigma|beta|period|shape",
                              "LP" = "latent_period",
                              "RP" = "removal_period",
                              "DP" = "detection_period",
                              "cov" = "cov_|^r_",
                              "fes" = "trial|donor|txd",
                              "weight" = "weight",
                              "all" = ".",
                              x),
                       negate = TRUE)
    })
    
    priors2 <- priors[parameter %in% pp]
    
    # Note: don't overwrite `use` or `type`. We want to patch the parameter
    # values, but the model determines whether they get used or not.
    params2$priors[parameter %in% priors2$parameter,
                   true_val := unlist(patch_vals)[parameter]]
    
    # Ensure val1 <= true_val <= val2
    params2$priors[, `:=`(val1 = pmin(val1, true_val),
                          val2 = pmax(val2, true_val))]
    
    
    # Patch from Trace file
    
    pars_patched <- "priors"
    
    # Covariance Matrices ----
    
    # We overwrite params2's Sigma_G and Sigma_E with new values, remembering
    # that r_*_** are correlations, not covariances, so they need to be
    # transformed first.
    if ("cov" %notin% skip_patches) {
        priors1 <- params2$priors$true_val |>
            setNames(params2$priors$parameter) |>
            as.list()
        patched_covariances <- FALSE
        
        params2$Sigma_G[] <- 0
        params2$Sigma_E[] <- 0
        params2$cov_G[] <- 0
        params2$cov_E[] <- 0
        
        out <- make_matrices_from_priors(priors1)
        
        # Copy the used values, discarding the rest (even if non-zero)
        used <- model_traits[str_chars(use_traits)]
        
        params2$Sigma_G[used, used] <- out$Sigma_G[used, used]
        params2$Sigma_E[used, used] <- out$Sigma_E[used, used]
        params2$cov_G[used, used] <- out$cov_G[used, used]
        params2$cov_E[used, used] <- out$cov_E[used, used]
        
        pars_patched <- c(pars_patched, "covariances")
    }
    
    params2$vars <- diag(params2$Sigma_G)
    params2$cors <- params2$Sigma_G[lower.tri(params2$Sigma_G)]

    
    # Rates and LPs ----
    
    if ("beta" %in% pp) {
        params2$r_beta <- patch_vals$beta
        pars_patched <- c(pars_patched, "beta")
    }
    if ("latent_period" %in% pp) {
        params2$latent_period <- patch_vals$latent_period
        pars_patched <- c(pars_patched, "latent period")
    }
    if ("LP_shape" %in% pp) {
        params2$LP_shape <- patch_vals$LP_shape
        params2$LP_scale <- patch_vals$latent_period / params2$LP_shape
        pars_patched <- c(pars_patched, "LP scale", "LP shape")
    }
    if ("detection_period" %in% pp) {
        params2$detection_period <- 1 / patch_vals$detection_period
        pars_patched <- c(pars_patched, "detection period")
    }
    if ("DP_shape" %in% pp) {
        params2$DP_shape <- patch_vals$DP_shape
        params2$DP_scale <- patch_vals$detection_period / params2$DP_shape
        pars_patched <- c(pars_patched, "DP scale", "DP shape")
    }
    if ("removal_period" %in% pp) {
        params2$removal_period <- patch_vals$removal_period
        pars_patched <- c(pars_patched, "removal period")
    }
    if ("RP_shape" %in% pp) {
        params2$RP_shape <- patch_vals$RP_shape
        params2$RP_scale <- patch_vals$removal_period / params2$RP_shape
        pars_patched <- c(pars_patched, "RP scale", "RP shape")
    }
    if ("sigma" %in% pp) {
        # params2$group_effect <- patch_vals$sigma
        pars_patched <- c(pars_patched, "group effect")
    }
    
    # Fixed effects ----
    
    # Handle nested Weight FEs
    weight_is_nested <- any(str_detect(names(patch_vals), "weight1"))
    weight <- if (weight_is_nested) c("weight1", "weight2") else "weight"

    fe_types <- c(
        if ("fes" %notin% skip_patches) c("trial", "donor", "txd"),
        if ("weight" %notin% skip_patches) weight
    )
    
    if (weight_is_nested) {
        params2$priors[str_starts(parameter, "weight_"), use := FALSE]
    }
    
    for (fe_type in fe_types) {
        for (fe_trait in model_traits) {
            # build name e.g. trial_l
            fe_name <- str_c(fe_type, "_", str_sub(fe_trait, 1, 1))
            if (fe_name %in% pp) {
                params2$fe_vals[fe_type, fe_trait] <- patch_vals[[fe_name]]
                pars_patched <- c(pars_patched, fe_name)
            }
        }
    }
    
    # BICI fixes
    
    if (params$sim_new_data %in% c("bici")) {
        params2 <- remove_bici_fes(params2)
    }
    
    if (resfile_exists) {
        if (FALSE) {
            # FIXME: this seems to be breaking, so just skip it
            
            # Ensure FEs match
            for (fe in fe_types) {
                fe_type <- str_c(fe, "_fe")
                if (fe_type %in% names(tmp_params)) {
                    params2[[fe_type]] <- tmp_params[[fe_type]]
                }
            }
            pars_patched <- c(pars_patched, "fixed effects")
        }
        
        # shape parameters
        # params2$LP_shape <- tmp_params$LP_shape
        # params2$DP_shape <- tmp_params$DP_shape
        # params2$RP_shape <- tmp_params$RP_shape
        
        # I think this is necessary?
        params2$time_step <- tmp_params$time_step
        params2$censor <- tmp_params$censor
        pars_patched <- c(pars_patched, "time_step", "censor")
    }
    
    # Ensure priors are wide enough to include true_vals
    widen_priors(params2$priors)
    
    if (msgs) message("- ", str_flatten_comma(pars_patched))
    
    params2
}
