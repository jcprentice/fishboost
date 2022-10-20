# This is a slower version but has the advantage of using very little memory
make_traits_from_pedigree2 <- function(pedigree, params) {
    message("Making trait values from pedigree ...")

    # extract parameters
    ntraits <- params$ntraits
    traitnames <- params$traitnames
    Sigma_G <- params$Sigma_G
    Sigma_E <- params$Sigma_E

    # drop traits where covariance = 0 (keep if small but non-zero)
    keep_traits <- traitnames[diag(Sigma_G) > 0]
    traitnames <- traitnames[keep_traits]
    ntraits <- length(keep_traits)
    Sigma_G <- Sigma_G[keep_traits, keep_traits]
    Sigma_E <- Sigma_E[keep_traits, keep_traits]


    # safer to get this directly from the pedigree
    ntotal <- nrow(pedigree)

    # create names of trait variants
    traitnames_BV <- paste0(traitnames, "_BV")

    traits <- copy(pedigree)

    # split into parents and progeny
    parent_ids  <- traits[sdp != "progeny", id]
    progeny_ids <- traits[sdp == "progeny", id]

    # set the genetic values (BVs) for the parents first
    traits[parent_ids, (traitnames_BV) := as.data.table(mvrnorm(.N, rep(0, ntraits), Sigma_G))]

    # set the genetic values of traits for the progeny
    for (id in progeny_ids) {
        # get ids of parents
        parents <- traits[id, c(sire, dam)]

        # get mean of parents' BVs
        parent_BVs <- traits[parents, lapply(.SD, mean), .SDcols = traitnames_BV]

        # generate new mvnorm values
        progeny_BV <- mvrnorm(1L, as.numeric(parent_BVs), Sigma_G / 2)

        # assign values to progeny
        set(traits, id, traitnames_BV, as.list(progeny_BV))
    }

    # generate environmental values
    traits_EV <- mvrnorm(ntotal, rep(0, ntraits), Sigma_E)

    # calculate traits as sum of genetic and environmental components
    # get log-normally distd values
    traits[, (traitnames)] <- exp(traits[, ..traitnames_BV] + traits_EV)

    traits
}
