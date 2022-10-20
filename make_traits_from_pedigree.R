# This is a faster version using arrays and no loops
make_traits_from_pedigree <- function(pedigree, params) {
    message("Making trait values from pedigree ...")

    # extract parameters
    all_traitnames <- params$all_traitnames
    traitnames     <- params$traitnames
    Sigma_G <- params$Sigma_G
    Sigma_E <- params$Sigma_E
    ntraits <- params$ntraits

    missing_traits <- setdiff(all_traitnames, traitnames)

    # skip everything if ntraits == 0 (otherwise things break)
    if (ntraits == 0) {
        traits <- copy(pedigree)
        traits[, (missing_traits) := 1]
        return(traits)
    }

    # safer to get these values directly from pedigree
    nparents <- pedigree[sdp != "progeny", .N]
    nprogeny <- pedigree[sdp == "progeny", .N]
    ntotal <- nrow(pedigree)

    # create names of trait variants
    traitnames_BV <- paste0(traitnames, "_BV")
    traitnames_EV <- paste0(traitnames, "_EV")

    # set the genetic values (BVs) for the parents first, and fix names
    parent_BVs  <- mvrnorm(nparents, rep(0, ntraits), Sigma_G)

    # get the parent IDs for the progeny in a convenient form
    sire_dam <- as.matrix(pedigree[, .(sire, dam)])

    # get mean parent BV
    progeny_ids <- pedigree[sdp == "progeny", id]

    sire_vals <- parent_BVs[sire_dam[progeny_ids, 1], ]
    dam_vals  <- parent_BVs[sire_dam[progeny_ids, 2], ]
    parent_means <- (sire_vals + dam_vals) / 2

    # generate progeny BVs ~ MvNorm(nprogeny, parent_means, Sigma_G / 2)
    progeny_BVs <- mvrnorm(nprogeny, rep(0, ntraits), Sigma_G / 2) + parent_means
    
    # Old method, may be more brittle depending on Sigma_G
    # L <- chol(Sigma_G / 2)
    # z <- matrix(rnorm(ntraits * nprogeny), nprogeny, ntraits)
    # progeny_BVs <- (z %*% L) + parent_means
    
    # generate environmental values
    EVs <- mvrnorm(ntotal, rep(0, ntraits), Sigma_E)
    dimnames(EVs)[[2]] <- traitnames_EV

    # generate log-normally dist'd phenotypes P ~ G + E
    BVs <- rbind(parent_BVs, progeny_BVs)
    dimnames(BVs)[[2]] <- traitnames_BV
    phenotypes <- exp(BVs + EVs)
    dimnames(phenotypes)[[2]] <- traitnames

    # combine pedigree, BVs, and phenotypes, to get traits
    traits <- cbind(pedigree, BVs, EVs, phenotypes)

    # we need these so the models run, but they don't do anything
    if (length(missing_traits) > 0) {
        traits[, (missing_traits) := 1]
    }

    # parents don't need phenotypes
    traits[sdp != "progeny", (all_traitnames) := NA]

    traits
}
