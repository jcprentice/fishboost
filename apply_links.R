# Take params and ensure that linked traits are correct by setting Sigma_G
apply_links <- function(params) {
    {
        all_traits       <- params$all_traits
        model_traits     <- params$model_traits
        use_traits       <- params$use_traits
        sim_trial_fe     <- params$sim_trial_fe
        sim_donor_fe     <- params$sim_donor_fe
        sim_txd_fe       <- params$sim_txd_fe
        sim_weight_fe    <- params$sim_weight_fe
        trial_fe         <- params$trial_fe
        donor_fe         <- params$donor_fe
        txd_fe           <- params$txd_fe
        weight_fe        <- params$weight_fe
        weight_is_nested <- params$weight_is_nested
        sim_link_traits  <- params$sim_link_traits
        sim_link_trial   <- params$sim_link_trial
        sim_link_donor   <- params$sim_link_donor
        sim_link_txd     <- params$sim_link_txd
        sim_link_weight  <- params$sim_link_weight
        sim_link_shapes  <- params$sim_link_shapes
        link_traits      <- params$link_traits
        link_trial       <- params$link_trial
        link_donor       <- params$link_donor
        link_txd         <- params$link_txd
        link_weight      <- params$link_weight
        link_shapes      <- params$link_shapes
        fe_vals          <- params$fe_vals
        Sigma_G          <- params$Sigma_G
        Sigma_E          <- params$Sigma_E
        cov_G            <- params$cov_G
        cov_E            <- params$cov_E
        LP_shape         <- params$LP_shape
        DP_shape         <- params$DP_shape
        RP_shape         <- params$RP_shape
    }

    # Need to make explicit copies because datatables work by reference
    params2 <- copy(params)

    # All changes made to priors also happen in params2
    priors  <- params2$priors

    # Grab the first letter of each model traitname for convenience
    sildt <- c("s", "i", "l", "d", "t")
    ntraits <- 5L


    # Individual Effects ----
    if (use_traits %notin% c("", "none")) {
        # Turn sim_link_traits into string vector
        new_idxs <- match(str_split_1(sim_link_traits, ""), sildt)

        # Copy rows and columns in Sigma and cov matrix
        params2$Sigma_G[] <- Sigma_G[new_idxs, new_idxs, drop = FALSE]
        params2$Sigma_E[] <- Sigma_E[new_idxs, new_idxs, drop = FALSE]
        params2$cov_G[] <- cov_G[new_idxs, new_idxs, drop = FALSE]
        params2$cov_E[] <- cov_E[new_idxs, new_idxs, drop = FALSE]

        # Be careful here in case any traits aren't used

        pars <- expand.grid("r_G_", sildt, sildt) |>
            apply(1, str_flatten) |>
            matrix(ntraits, ntraits) |>
            {\(x) {
                diag(x) <- str_replace(diag(x), "r", "cov")
                x
            }}() |>
            {\(x) {
                y <- c(diag(x), t(x)[lower.tri(x)])
                c(y, str_replace(y, "G", "E"))
            }}()

        parvals <- data.table(par = pars,
                              val = c(diag(Sigma_G), Sigma_G[lower.tri(Sigma_G)],
                                      diag(Sigma_E), Sigma_E[lower.tri(Sigma_E)]))

        priors[parameter %in% parvals$par,
               true_val := fifelse(use, parvals$val, true_val)]
    }


    # Fixed Effects ----
    if (sim_trial_fe %notin% c("", "none")) {
        # Copy FE values across, overwriting previous value
        zt <- match(str_split_1(sim_link_trial, ""), sildt)
        fe_vals["trial", ] <- fe_vals["trial", zt]
    }

    if (trial_fe %notin% c("", "none")) {
        # sieve "SILDR" for trial_fe, map, and convert into list of priors
        trial_priors <- str_c("trial_",
                              str_split_1(link_trial, "") |>
                                  intersect(sildt) |>
                                  intersect(str_split_1(trial_fe, "")))
        priors[str_starts(parameter, "trial_"), use := parameter %in% trial_priors]
    }

    if (sim_donor_fe %notin% c("", "none")) {
        zd <- match(str_split_1(sim_link_donor, ""), sildt)
        fe_vals["donor", ] <- fe_vals["donor", zd]
    }

    if (donor_fe %notin% c("", "none")) {
        donor_priors <- str_c("donor_",
                              str_split_1(link_donor, "") |>
                                  intersect(sildt) |>
                                  intersect(str_split_1(donor_fe, "")))
        priors[str_starts(parameter, "donor_"), use := parameter %in% donor_priors]
    }

    if (sim_txd_fe %notin% c("", "none")) {
        zd <- match(str_split_1(sim_link_txd, ""), sildt)
        fe_vals["txd", ] <- fe_vals["txd", zd]
    }

    if (txd_fe %notin% c("", "none")) {
        txd_priors <- str_c("txd_",
                            str_split_1(link_txd, "") |>
                                intersect(sildt) |>
                                intersect(str_split_1(txd_fe, "")))
        priors[str_starts(parameter, "txd_"), use := parameter %in% txd_priors]
    }

    if (sim_weight_fe %notin% c("", "none")) {
        zd <- match(str_split_1(sim_link_weight, ""), sildt)
        wts <- str_subset(rownames(fe_vals), "weight")
        fe_vals[wts, ] <- fe_vals[wts, zd]
    }

    if (weight_fe %notin% c("", "none")) {
        link_weight_v <- str_split_1(link_weight, "") |>
            intersect(sildt) |>
            intersect(str_split_1(weight_fe, ""))
        if (weight_is_nested) {
            weight1_priors <- str_c("weight1_", link_weight_v)
            weight2_priors <- str_c("weight2_", link_weight_v)
            priors[str_starts(parameter, "weight1_"), use := parameter %in% weight1_priors]
            priors[str_starts(parameter, "weight2_"), use := parameter %in% weight2_priors]
        } else {
            weight_priors <- str_c("weight_", link_weight_v)
            priors[str_starts(parameter, "weight_"), use := parameter %in% weight_priors]
        }
    }

    # Copy updated values back into params2
    params2$fe_vals <- fe_vals

    # Shape parameters ----
    shapes <- c(l = LP_shape, d = DP_shape, t = RP_shape)
    shape <- shapes[str_split_1(sim_link_shapes, "")]
    params2[c("LP_shape", "DP_shape", "RP_shape")] <- shape

    # Set values in priors (be careful not to accidentally turn them back on if the type is "exp")
    priors[parameter == "LP_shape", use := use && str_detect(link_shapes, "l")]
    priors[parameter == "DP_shape", use := use && str_detect(link_shapes, "d")]
    priors[parameter == "RP_shape", use := use && str_detect(link_shapes, "t")]

    params2
}

