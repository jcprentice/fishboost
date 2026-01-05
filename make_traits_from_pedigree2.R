# This is a slower version but has the advantage of using very little memory
make_traits_from_pedigree2 <- function(popn, params) {
    message("Making trait values from pedigree ...")

    # extract parameters
    n_traits <- params$n_traits
    traitnames <- params$traitnames
    cov_G <- params$cov_G
    cov_E <- params$cov_E

    # drop popn where covariance = 0 (keep if small but non-zero)
    keep_traits <- traitnames[diag(cov_G) > 0]
    traitnames <- traitnames[keep_traits]
    n_traits <- length(keep_traits)
    cov_G <- cov_G[keep_traits, keep_traits]
    cov_E <- cov_E[keep_traits, keep_traits]


    # safer to get this directly from the pedigree
    ntotal <- nrow(popn)

    # create names of trait variants
    traitnames_GV <- str_c(traitnames, "_g")

    popn2 <- copy(popn)

    # split into parents and progeny
    parent_ids  <- popn2[sdp != "progeny", id]
    progeny_ids <- popn2[sdp == "progeny", id]

    # set the genetic values (BVs) for the parents first
    popn2[parent_ids, (traitnames_GV) := as.data.table(mvrnorm(.N, rep(0, n_traits), cov_G))]

    # set the genetic values of popn for the progeny
    walk(progeny_ids, \(id) {
        # get ids of parents
        parents <- popn2[id, c(sire, dam)]

        # get mean of parents' BVs
        parent_GVs <- popn2[parents, map(.SD, mean), .SDcols = traitnames_GV]

        # generate new mvnorm values
        progeny_GV <- mvrnorm(1L, as.numeric(parent_GVs), cov_G / 2)

        # assign values to progeny
        set(popn2, id, traitnames_GV, as.list(progeny_GV))
    }

    # generate environmental values
    popn_EV <- mvrnorm(ntotal, rep(0, n_traits), cov_E)

    # calculate popn as sum of genetic and environmental components
    # get log-normally distd values
    popn2[, (traitnames)] <- exp(popn2[, ..traitnames_GV] + popn_EV)

    # Add weights
    if ("weight" %notin% names(popn)) {
        popn2[, weight := 1]
        popn2[trial == 1, weight := rlnorm(.N, 3.430757, 0.2899917)]
        popn2[trial == 2, weight := rlnorm(.N, 4.457181, 0.3330279)]
    }

    popn2
}
