# Further details of this sampling procedure can be seen here:
# https://en.wikipedia.org/wiki/Multivariate_normal_distribution#Drawing_values_from_the_distribution

make_traits_from_grm <- function(GRM, pedigree, params) {
    message("Making trait values from GRM ...")

    ntotal     <- nrow(pedigree)
    ntraits    <- params$ntraits
    Sigma_G    <- params$Sigma_G
    Sigma_E    <- params$Sigma_E
    traitnames <- params$traitnames

    traitnames_BV <- paste0(traitnames, "_BV")

    # Generate genetic values (BVs)
    KGG <- kronecker(Matrix(GRM), Sigma_G)

    # note that R's Cholesky fn gives a U matrix, whereas the wikipedia
    # algorithm wants an L matrix, so need to take the transpose
    L <- t(chol(KGG))

    z <- rnorm(ntotal * ntraits)
    traits_BV <- matrix(L %*% z,
                        nrow = ntotal, ncol = ntraits,
                        byrow = TRUE,
                        dimnames = list(NULL, traitnames_BV))

    # Generate environmental values
    # much easier as these don't need to account for relatedness
    traits_EV <- mvrnorm(ntotal, numeric(ntraits), Sigma_E)

    # traits starts as a copy of pedigree
    traits <- copy(pedigree)

    # trait = exp(genetic + environment components)
    traits[, (traitnames_BV)] <- as.data.table(traits_BV)
    traits[, (traitnames)] <- as.data.table(exp(traits_BV + traits_EV))

    traits
}
