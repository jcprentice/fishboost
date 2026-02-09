library(data.table)
library(stringr)

#' Apply trial, donor, and weight FEs to popn.
#'
#' @description `apply_fixed_effects()` takes a population data.table and a
#'   parameters list and applies all the fixed effects defined by `sim_X_fe`
#'   according to the FE values contained in `fe_vals`.
#'
#' @param popn A population data.table
#' @param params A parameters list
#'
#' @returns A new population file with FEs applied to phenotypes

apply_fixed_effects <- function(popn, params) {
    {
        model_traits     <- params$model_traits
        sim_trial_fe     <- params$sim_trial_fe
        sim_donor_fe     <- params$sim_donor_fe
        sim_txd_fe       <- params$sim_txd_fe
        sim_weight_fe    <- params$sim_weight_fe
        fe_vals          <- params$fe_vals
        use_weight       <- params$use_weight
        weight_is_nested <- params$weight_is_nested
    }

    message("Applying fixed effects to popn...")

    popn2 <- copy(popn)

    # Some traits don't have GVs and EVs, so we give placeholders with value 0

    traits_g <- str_c(model_traits, "_g")
    traits_e <- str_c(model_traits, "_e")
    missing_GEVs <- setdiff(c(traits_g, traits_e),
                            names(popn2))
    popn2[, (missing_GEVs) := 0]
    
    # Choose how to recentre the weights, log (default) or linear
    recentre_f <- get(if (use_weight == "log") "log_recentre" else "recentre")
    
    N1 <- popn2[trial == 1, .N]
    N2 <- popn2[trial == 2, .N]
    
    mat <- popn2[sdp == "progeny",
                 .(trial = recentre(trial == 2),
                   donor = recentre(donor == 1),
                   txd   = recentre(trial == 2 & donor == 1),
                   weight  = recentre_f(weight),
                   weight1 = c(recentre_f(weight[trial == 1]), rep(0, N2)),
                   weight2 = c(rep(0, N1), recentre_f(weight[trial == 2]))
                 )] |>
        as.matrix()
    
    # Zero out all fe_vals that aren't in sim_X_fe
    fe_vals["trial",   sildt %notin% str_chars(sim_trial_fe)]  <- 0
    fe_vals["donor",   sildt %notin% str_chars(sim_donor_fe)]  <- 0
    fe_vals["txd",     sildt %notin% str_chars(sim_txd_fe)]    <- 0
    fe_vals["weight",  sildt %notin% str_chars(sim_weight_fe)] <- 0
    fe_vals["weight1", sildt %notin% str_chars(sim_weight_fe)] <- 0
    fe_vals["weight2", sildt %notin% str_chars(sim_weight_fe)] <- 0
    
    if (weight_is_nested) {
        fe_vals["weight", ] <- 0
    } else {
        fe_vals[c("weight1", "weight2"), ] <- 0
    }
    
    
    mat_fe <- mat %*% fe_vals
    
    mat_g <- popn2[sdp == "progeny", ..traits_g] |> as.matrix()
    mat_e <- popn2[sdp == "progeny", ..traits_e] |> as.matrix()
    
    # Calculate the phenotype
    mat_pt <- exp(mat_g + mat_e + mat_fe)
    colnames(mat_pt) <- model_traits
    
    popn2[sdp == "progeny", (model_traits) := as.data.frame(mat_pt)]

    popn2[, (missing_GEVs) := NULL]

    popn2
}
