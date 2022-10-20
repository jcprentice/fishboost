# Spearman's rank correlation
spearmans_rc <- function(x, y) {
    rx <- rank(x)
    ry <- rank(y)
    n <- length(x)
    1 - 6 * sum((rx - ry)^2) / (n * (n^2 - 1))
}


# Return EBVs ranked, show rank distance vs true ranks
get_ranks <- function(pop, estimated_BVs, params, verbose = FALSE) {
    {
        traitnames <- params$traitnames
        sire_version <- params$sire_version
    }

    # make a copy of this so we can add missing columns
    pop2 <- copy(pop)

    sir <- c("susceptibility", "infectivity", "recoverability")
    used <- sir %in% traitnames

    # individual effect names (what SIRE 2.1 calls them)
    
    if (sire_version == "2.2") {
        ie_names <- c("s_a EBV", "i_a EBV", "r_a EBV")
    } else {
        ie_names <- c("s_a", "i_a", "r_a")
    }

    setnames(estimated_BVs, ie_names, sir, skip_absent = TRUE)

    # add missing columns
    if (!used[1]) {
        estimated_BVs[, susceptibility := 0]
        pop2[, susceptibility_BV := 0]
    }
    if (!used[2]) {
        estimated_BVs[, infectivity := 0]
        pop2[, infectivity_BV := 0]
    }
    if (!used[3]) {
        estimated_BVs[, recoverability := 0]
        pop2[, recoverability_BV := 0]
    }

    ids <- seq.int(params$nsires)

    BVs <- estimated_BVs[ids, .(id,
                                est_sus = susceptibility,
                                est_inf = infectivity,
                                est_rec = recoverability)]

    # Can't make comparisons if using Fishboost data
    if (params$use_fb_data == FALSE) {
        BVs2 <- pop2[ids, .(id,
                           sus = susceptibility_BV,
                           inf = infectivity_BV,
                           rec = recoverability_BV)]

        BVs <- merge(BVs2, BVs, by = "id")

        if (verbose) {
            message("Spearman's rank distance for:")
            message(" - susceptibility = ", signif(cor.test(BVs$sus, BVs$est_sus, method = "spearman")$estimate, 3))
            message(" - infectivity    = ", signif(cor.test(BVs$inf, BVs$est_inf, method = "spearman")$estimate, 3))
            message(" - recoverability = ", signif(cor.test(BVs$rec, BVs$est_rec, method = "spearman")$estimate, 3))
        }
    }

    ranks <- BVs[, lapply(.SD, rank)]
}

