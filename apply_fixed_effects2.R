library(data.table)
library(stringr)

# Apply trial and donor FEs, by making a copy of popn, copying the GV,
# applying the FE to the relevant individuals, normalising the population, and
# recalculating the trait

apply_fixed_effects2 <- function(popn, params) {
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
    missing_GEVs <- setdiff(all_GEVs, names(popn))
    popn2[, (missing_GEVs) := 0]

    # the set of donors and trial2 individuals
    popn2[, `:=`(trial2 = trial == 2,
                 donor2 = donor == 1,
                 txd2 = trial == 2 & donor == 1)]

    # looking for either log_recentre() or recentre()
    recentre_f <- get(if (use_weight == "log") "log_recentre" else "recentre")
    
    # TODO: remake this as a matrix operation if possible
    
    walk(model_traits, \(trait) {
        # trait <- model_traits[[1]]
        trait1 <- str_1st(trait)
        
        # str_detect() returns TRUE or FALSE, which is effectively 1 or 0 when
        # multiplied by a FE value, so an unused trait becomes 0.
        trial_val   <- str_detect(sim_trial_fe,  trait1) * fe_vals["trial",   trait]
        donor_val   <- str_detect(sim_donor_fe,  trait1) * fe_vals["donor",   trait]
        txd_val     <- str_detect(sim_txd_fe,    trait1) * fe_vals["txd",     trait]
        weight_val  <- str_detect(sim_weight_fe, trait1) * fe_vals["weight",  trait]
        weight1_val <- str_detect(sim_weight_fe, trait1) * fe_vals["weight1", trait]
        weight2_val <- str_detect(sim_weight_fe, trait1) * fe_vals["weight2", trait]

        popn2[, tmp := trial2 * trial_val + donor2 * donor_val + txd2 * txd_val]

        if (DEBUG && msgs) {
            message(trait)
            ls() |> str_subset("_val$") |> mget() |> unlist() |> print() |> capture_message()
        }

        # Weight needs a bit more care
        if (weight_is_nested) {
            popn2[trial == 1, tmp := tmp + weight1_val * recentre_f(weight)]
            popn2[trial == 2, tmp := tmp + weight2_val * recentre_f(weight)]
        } else {
            popn2[trial %in% 1:2, tmp := tmp + weight_val * recentre_f(weight)]
        }

        popn2[, tmp := tmp - mean(tmp, na.rm = TRUE)]

        popn2[, (trait) := {
            gv <- get(str_c(trait, "_g"))
            ev <- get(str_c(trait, "_e"))
            exp(gv + ev + tmp)
        }]
    })

    # Clean up
    popn2[, c("trial2", "donor2", "txd2", "tmp") := NULL]
    popn2[, (missing_GEVs) := NULL]

    popn2
}
