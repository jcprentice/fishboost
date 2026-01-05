{
    library(data.table)
    library(purrr)
}

# If the protocol touches any of the TPs or shape parameters, make sure that
# these changes are reflected in the rate parameters and the priors.

tidy_up_periods <- function(params, protocol = NULL) {
    params2      <- copy(params)
    sim_new_data <- params$sim_new_data
    
    TPs <- if (is.null(protocol)) {
        c("latent_period",    "LP_shape",
          "detection_period", "DP_shape",
          "removal_period",   "RP_shape")
    } else {
        names(protocol)
    }
    
    # Fix TP scales
    params2$LP_scale <- with(params2, latent_period    / LP_shape)
    params2$DP_scale <- with(params2, detection_period / DP_shape)
    params2$RP_scale <- with(params2, removal_period   / RP_shape)
    
    # Update the TP priors
    walk(TPs,
         ~ params2$priors[parameter == .x, `:=`(
             val1 = min(params2[[.x]], val1),
             val2 = max(params2[[.x]], val2),
             true_val = params2[[.x]]
         )]
    )
    
    params2
}
