# Spearman's rank correlation
spearmans_rc <- function(x, y) {
    rx <- rank(x)
    ry <- rank(y)
    n <- length(x)
    1 - 6 * sum((rx - ry)^2) / (n * (n^2 - 1))
}


# Return EBVs ranked, show rank distance vs true ranks
get_ranks <- function(popn, estimated_BVs, params, verbose = FALSE) {
    {
        use_traits   <- params$use_traits
        sire_version <- params$sire_version
    }

    # make a copy of this so we can add missing columns
    popn2 <- copy(popn)

    # individual effect names (what SIRE 2.1 calls them)
    
    ie_names <- if (sire_version == "bici") {
        c("sg", "ig", "tg")
    } else {
        c("s_g EBV", "i_g EBV", "t_g EBV")
    }

    setnames(estimated_BVs, ie_names, c("sus", "inf", "tol"), skip_absent = TRUE)

    # add missing columns
    if (!str_detect(use_traits, "s")) {
        estimated_BVs[, sus := 0]
        popn2[, sus_g := 0]
    }
    if (!str_detect(use_traits, "i")) {
        estimated_BVs[, inf := 0]
        popn2[, inf_g := 0]
    }
    if (!str_detect(use_traits, "t")) {
        estimated_BVs[, tol := 0]
        popn2[, tol_g := 0]
    }

    ids <- seq_len(params$nsires)

    BVs <- estimated_BVs[ids, .(id,
                                est_sus = sus,
                                est_inf = inf,
                                est_tol = tol)]

    # Can't make comparisons if using Fishboost data
    if (params$sim_new_data != "no") {
        BVs2 <- popn2[ids, .(id,
                           sus = sus_g,
                           inf = inf_g,
                           tol = tol_g)]

        BVs <- merge(BVs2, BVs, by = "id")

        if (verbose) {
            message("Spearman's rank distance for:")
            message(" - susceptibility = ", signif(cor.test(BVs$sus, BVs$est_sus, method = "spearman")$estimate, 3))
            message(" - infectivity    = ", signif(cor.test(BVs$inf, BVs$est_inf, method = "spearman")$estimate, 3))
            message(" - tolerance      = ", signif(cor.test(BVs$tol, BVs$est_tol, method = "spearman")$estimate, 3))
        }
    }

    ranks <- BVs[, map(.SD, rank)]
}

