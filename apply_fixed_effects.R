# Apply trial and donor FEs, by making a copy of traits, copying the BV,
# applying the FE to the relevant individuals, normalising the population, and
# recalculating the trait

apply_fixed_effects <- function(traits, params) {
    {
        traitnames     <- params$traitnames
        all_traitnames <- params$all_traitnames
        sim_trial_fe   <- params$sim_trial_fe
        sim_donor_fe   <- params$sim_donor_fe
        sim_txd_fe     <- params$sim_txd_fe
        fe_vals        <- params$fe_vals
    }

    traits2 <- copy(traits)

    # Some traits don't have BVs and EVs, need to generate fake ones = 0
    missing <- setdiff(all_traitnames, traitnames)
    missing_Vs <- c(paste0(missing, "_BV"), paste0(missing, "_EV"))
    traits2[, (missing_Vs) := 0]

    # the set of donors and trial2 individuals
    traits2[, `:=`(trial2 = trial == 2,
                   donor2 = donor == 1,
                   txd2 = trial == 2 & donor == 1)]

    for (trait in all_traitnames) {
        trial_val <- grepl(substr(trait, 1, 1), sim_trial_fe) * fe_vals["trial", trait]
        donor_val <- grepl(substr(trait, 1, 1), sim_donor_fe) * fe_vals["donor", trait]
        txd_val   <- grepl(substr(trait, 1, 1), sim_txd_fe)   * fe_vals["txd",    trait]

        traits2[, tmp := trial2 * trial_val + donor2 * donor_val + txd2 * txd_val]
        traits2[, tmp := tmp - mean(tmp, na.rm = TRUE)]

        trait_BV <- paste0(trait, "_BV")
        trait_EV <- paste0(trait, "_EV")
        traits2[, (trait) := exp(get(trait_BV) + get(trait_EV) + tmp)]
    }

    # Clean up
    traits2[, c("trial2", "donor2", "txd2", "tmp") := NULL]
    traits2[, (missing_Vs) := NULL]

    traits2
}
