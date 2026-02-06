library(data.table)
library(stringr)

# Apply trial and donor FEs, by making a copy of popn, copying the GV,
# applying the FE to the relevant individuals, normalising the population, and
# recalculating the trait

apply_fixed_effects <- function(popn, params) {
    # message("Applying fixed effects ...")

    {
        model_traits     <- params$model_traits
        sim_trial_fe     <- params$sim_trial_fe
        sim_donor_fe     <- params$sim_donor_fe
        sim_txd_fe       <- params$sim_txd_fe
        sim_weight_fe    <- params$sim_weight_fe
        fe_vals          <- params$fe_vals
        use_weight       <- params$use_weight
        weight_is_nested <- params$weight_is_nested
        msgs             <- params$msgs
        DEBUG            <- params$DEBUG
    }

    popn2 <- copy(popn)

    # Some traits don't have GVs and EVs, need to generate fake ones = 0
    all_GEVs <- c(str_c(model_traits, "_g"), str_c(model_traits, "_e"))
    missing_GEVs <- setdiff(all_GEVs, names(popn2))
    popn2[, (missing_GEVs) := 0]
    
    # looking for either log_recentre() or recentre()
    recentre_f <- get(if (use_weight == "log") "log_recentre" else "recentre")
    
    mat <- popn2[sdp == "progeny",
                 .(trial = recentre(trial == 2),
                   donor = recentre(donor == 1),
                   txd   = recentre(trial == 2 & donor == 1),
                   weight  = recentre_f(weight),
                   weight1 = c(recentre_f(weight[trial == 1]), rep(0, sum(trial == 1))),
                   weight2 = c(rep(0, sum(trial == 2)), recentre_f(weight[trial == 2]))
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
    
    traits_g <- str_c(model_traits, "_g")
    traits_e <- str_c(model_traits, "_e")
    mat_g <- popn2[sdp == "progeny", ..traits_g] |> as.matrix()
    mat_e <- popn2[sdp == "progeny", ..traits_e] |> as.matrix()
    
    # Calculate the phenotype
    mat_pt <- exp(mat_g + mat_e + mat_fe)
    colnames(mat_pt) <- model_traits
    
    popn2[sdp == "progeny", (model_traits) := as.data.frame(mat_pt)]

    popn2[, (missing_GEVs) := NULL]

    popn2
}
