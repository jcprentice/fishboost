library(MASS)

# This is a faster version using arrays and no loops
make_traits_from_pedigree <- function(pedigree, params) {
    message("Making trait values from pedigree ...")

    # Extract parameters
    {
        model_traits <- params$model_traits
        use_traits   <- params$use_traits
        traitnames   <- model_traits[intersect(str_chars(use_traits),
                                               names(model_traits))]
        cov_G        <- params$cov_G[traitnames, traitnames]
        cov_E        <- params$cov_E[traitnames, traitnames]
        n_traits     <- str_length(use_traits)
        use_weight   <- params$use_weight
    }

    # Keep track of traits we haven't used, for those the genotypic value will
    # be set to 0 and the phenotypic value will be set to 1.
    missing_traits <- setdiff(model_traits, traitnames)

    # Skip everything if n_traits == 0 (otherwise things break)
    if (n_traits == 0) {
        popn <- copy(pedigree)
        missing_traits_GV <- str_c(missing_traits, "_g")
        missing_traits_EV <- str_c(missing_traits, "_e")
        # popn[, trial := 1L] # what is this doing here?
        popn[, (missing_traits_GV) := 0]
        popn[, (missing_traits_EV) := 0]
        popn[, (missing_traits) := 1]
        return(popn)
    }

    # Safer to get these values directly from pedigree
    nparents <- pedigree[sdp != "progeny", .N]
    nprogeny <- pedigree[sdp == "progeny", .N]
    ntotal <- nrow(pedigree)

    # Create names of trait variants
    traitnames_GV <- str_c(traitnames, "_g")
    traitnames_EV <- str_c(traitnames, "_e")
    
    # Check for positive definite matrix
    if (any(eigen(cov_G)$values <= 0)) {
        message("- Matrix is not positive definite, shrinking correlations")
        n <- 0
        while (n < 50 && any(eigen(cov_G)$values < 0)) {
            message("shrinking")
            n <- n + 1
            dcov <- diag(cov_G)
            cov_G <- 0.9 * cov_G
            diag(cov_G) <- dcov
        }
    }

    # Set the genetic values (GVs) for the parents first, and fix names
    parent_GVs <- mvrnorm(nparents, rep(0, n_traits), cov_G)

    # Get the parent IDs for the progeny in a convenient form
    sire_dam <- as.matrix(pedigree[, .(sire, dam)])

    # Get mean parent GV
    progeny_ids <- pedigree[sdp == "progeny", id]

    sire_vals <- parent_GVs[sire_dam[progeny_ids, 1], ]
    dam_vals  <- parent_GVs[sire_dam[progeny_ids, 2], ]
    
    # Correct for missing values in pedigree
    na_sire_ids <- pedigree[sdp == "progeny" & is.na(sire), id - nparents]
    if (length(na_sire_ids) > 0) {
        sire_vals[na_sire_ids, ] <- mvrnorm(length(na_sire_ids), rep(0, n_traits), cov_G)
    }
    
    na_dam_ids <- pedigree[sdp == "progeny" & is.na(dam), id - nparents]
    if (length(na_dam_ids) > 0) {
        dam_vals[na_dam_ids, ] <- mvrnorm(length(na_dam_ids), rep(0, n_traits), cov_G)
    }
    
    parent_means <- (sire_vals + dam_vals) / 2

    # Generate progeny GVs ~ MvNorm(nprogeny, parent_means, cov_G / 2)
    progeny_GVs <- mvrnorm(nprogeny, rep(0, n_traits), cov_G / 2) + parent_means
    
    # Old method, may be more brittle depending on cov_G
    # L <- chol(cov_G / 2)
    # z <- matrix(rnorm(n_traits * nprogeny), nprogeny, n_traits)
    # progeny_GVs <- (z %*% L) + parent_means
    
    # Generate environmental values
    EVs <- mvrnorm(ntotal, rep(0, n_traits), cov_E)
    colnames(EVs) <- traitnames_EV
    
    # Generate log-normally dist'd phenotypes P ~ G + E
    GVs <- rbind(parent_GVs, progeny_GVs)
    colnames(GVs) <- traitnames_GV
    
    # We need to shift the traits by the mean in order for the exp() values to
    # have mean 1.
    GVs <- GVs - apply(GVs, 2, var) / 2
    EVs <- EVs - apply(EVs, 2, var) / 2
    
    phenotypes <- exp(GVs + EVs)
    colnames(phenotypes) <- traitnames

    # Combine pedigree, GVs, and phenotypes, to get traits
    popn <- cbind(pedigree, GVs, EVs, phenotypes)

    # we need these so the models run, but they don't do anything
    if (length(missing_traits) > 0) {
        popn[, (missing_traits) := 1]
    }
    
    # parents don't need phenotypes
    popn[sdp != "progeny", (model_traits) := NA]

    popn
}
