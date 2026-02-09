#' Ensure linked traits are correctly applied
#'
#' @description Take `params` and ensure that linked traits are correct by
#'   setting `Sigma_G`, `priors`, and `fe_vals.`
#'
#' @param params The list of model parameters
#' @returns A new modified model parameters

apply_links <- function(params) {
    {
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
        LP_shape         <- params$LP_shape %||% "ldt"
        DP_shape         <- params$DP_shape %||% "ldt"
        RP_shape         <- params$RP_shape %||% "ldt"
    }

    # Make explicit copies because data.tables work by reference
    params2 <- copy(params)
    
    # All changes made to priors also happen in params2
    priors  <- params2$priors

    # Get the first letter of each model traitname for convenience
    sildt <- c("s", "i", "l", "d", "t")
    ntraits <- 5L

    get_idxs <- function(s) match(str_chars(s), sildt)

    # Individual Effects ----
    if (use_traits %notin% c("", "none")) {
        # Turn sim_link_traits into string vector
        new_idxs <- get_idxs(sim_link_traits)

        # Copy rows and columns in Sigma and cov matrix
        params2$Sigma_G[] <- Sigma_G[new_idxs, new_idxs, drop = FALSE]
        params2$Sigma_E[] <- Sigma_E[new_idxs, new_idxs, drop = FALSE]
        params2$cov_G[] <- cov_G[new_idxs, new_idxs, drop = FALSE]
        params2$cov_E[] <- cov_E[new_idxs, new_idxs, drop = FALSE]

        # Be careful here in case any traits aren't used
        set_diag_to_cov <- function(x) {
            diag(x) <- str_replace(diag(x), "r", "cov")
            x
        }
        get_LT_pars <- function(x) c(diag(x), t(x)[lower.tri(x)])
        duplicate_GE <- function(x) c(x, str_replace(x, "G", "E"))
        
        pars <- expand.grid("r_G_", sildt, sildt) |>
            apply(1, str_flatten) |>
            matrix(ntraits, ntraits) |>
            set_diag_to_cov() |>
            get_LT_pars() |>
            duplicate_GE()

        parvals <- data.table(par = pars,
                              val = c(diag(Sigma_G),
                                      Sigma_G[lower.tri(Sigma_G)],
                                      diag(Sigma_E),
                                      Sigma_E[lower.tri(Sigma_E)]))

        priors[parameter %in% parvals$par,
               true_val := fifelse(use, parvals$val, true_val)]
    }

    get_priors <- function(fe = "trial") {
        str  <- get(str_c(fe, "_fe"))   |> str_chars() # ildtt
        link <- get(str_c("link_", fe)) |> str_chars() # sittt
        str_c(fe, "_", intersect(str, link)) # trial_i, trial_t
    }

    # Fixed Effects ----
    if (sim_trial_fe %notin% c("", "none")) {
        # Copy FE values across, overwriting previous value
        fe_vals["trial", ] <- fe_vals["trial", get_idxs(sim_link_trial)]
    }

    if (trial_fe %notin% c("", "none")) {
        # sieve "sildt" for trial_fe, map, and convert into list of priors
        priors[str_starts(parameter, "trial_"),
               use := parameter %in% get_priors("trial")]
    }

    if (sim_donor_fe %notin% c("", "none")) {
        fe_vals["donor", ] <- fe_vals["donor", get_idxs(sim_link_donor)]
    }

    if (donor_fe %notin% c("", "none")) {
        priors[str_starts(parameter, "donor_"),
               use := parameter %in% get_priors("donor")]
    }

    if (sim_txd_fe %notin% c("", "none")) {
        fe_vals["txd", ] <- fe_vals["txd", get_idxs(sim_link_txd)]
    }

    if (txd_fe %notin% c("", "none")) {
        priors[str_starts(parameter, "txd_"),
               use := parameter %in% get_priors("txd")]
    }

    if (sim_weight_fe %notin% c("", "none")) {
        wts <- str_subset(rownames(fe_vals), "weight")
        fe_vals[wts, ] <- fe_vals[wts, get_idxs(sim_link_weight)]
    }
    
    if (weight_fe %notin% c("", "none")) {
        wt <- get_priors("weight")
        if (weight_is_nested) {
            priors[str_starts(parameter, "weight1_"),
                   use := parameter %in% str_replace(wt, "_", "1_")]
            priors[str_starts(parameter, "weight2_"),
                   use := parameter %in% str_replace(wt, "_", "2_")]
        } else {
            priors[str_starts(parameter, "weight_"), use := parameter %in% wt]
        }
    }

    # Copy updated values back into params2
    params2$fe_vals <- fe_vals

    # Shape parameters ----
    shapes <- c(l = LP_shape, d = DP_shape, t = RP_shape)
    shape <- shapes[str_chars(sim_link_shapes)]
    params2[c("LP_shape", "DP_shape", "RP_shape")] <- shape

    # Set values in priors (be careful not to accidentally turn them back on if the type is "exp")
    priors[parameter == "LP_shape", use := use && str_detect(link_shapes, "l")]
    priors[parameter == "DP_shape", use := use && str_detect(link_shapes, "d")]
    priors[parameter == "RP_shape", use := use && str_detect(link_shapes, "t")]

    params2
}

