# Further details of this sampling procedure can be seen here:
# https://en.wikipedia.org/wiki/Multivariate_normal_distribution#Drawing_values_from_the_distribution

# Also note that the same method can be used with the precision matrix, Sigma^-1
# https://math.stackexchange.com/questions/1845132/sample-points-from-a-multivariate-normal-distribution-using-only-the-precision-m

make_traits_from_grm <- function(popn, params) {
    message("Making trait values from GRM ...")
    
    {
        cov_G        <- params$cov_G
        cov_E        <- params$cov_E
        use_traits   <- params$use_traits
        model_traits <- params$model_traits
        use_weight   <- params$use_weight
        use_grm      <- params$use_grm
        setup        <- params$setup
    }
    
    message(str_glue(" - using '{use_grm}'"))
    
    GRM <- if (str_starts(use_grm, "A")) {
        make_grm(popn, "A")
    } else if (str_starts(use_grm, "H")) {
        # just want Hx_inv, not _nz
        fread(str_glue("fb_data/{x}_inv_{y}.tsv",
                       x = str_split_i(use_grm, "_", 1), 
                       y = str_remove(setup, "fb_")),
              header = TRUE) |>
            as.matrix() |>
            solve()
    } else {
        stop(" - Cannot find GRM!")
    }
    
    # This should stop a failure at chol(KGG) if a diagonal element is 0
    # diag(cov_G) <- pmax(diag(cov_G), 1e-6)
    # diag(cov_E) <- pmax(diag(cov_E), 1e-6)
    
    traits <- model_traits[str_detect(use_traits, names(model_traits))]
    n_used_traits <- length(traits)
    unused_traits <- setdiff(model_traits, traits)
    traitnames_GV <- str_c(traits, "_g")
    traitnames_EV <- str_c(traits, "_e")
    
    
    # Generate genetic values (GVs)
    KGG <- kronecker(GRM, cov_G[traits, traits])
    
    # Note that R's Cholesky fn gives a U matrix, whereas the wikipedia
    # algorithm wants an L matrix
    L <- chol(KGG)
    
    ntotal <- nrow(popn)
    z <- rnorm(ntotal * n_used_traits)
    traits_GV <- matrix(z %*% L,
                        nrow = ntotal, ncol = n_used_traits,
                        byrow = TRUE,
                        dimnames = list(NULL, traitnames_GV))
    
    # Generate environmental values
    # much easier as these don't need to account for relatedness
    traits_EV <- mvrnorm(ntotal, rep(0, n_used_traits), cov_E[traits, traits])
    
    # don't destroy original popn
    popn2 <- copy(popn)
    
    # trait = exp(genetic + environment components)
    popn2[, (traitnames_GV)] <- as.data.table(traits_GV)
    popn2[, (traitnames_EV)] <- as.data.table(traits_EV)
    popn2[, (traits)] <- as.data.table(exp(traits_GV + traits_EV))
    
    # we need these so the models run, but they don't do anything
    if (length(unused_traits) > 0) {
        popn2[, (unused_traits) := 1]
    }
    
    # parents don't need phenotypes
    popn2[sdp != "progeny", (model_traits) := NA]
    
    # Add weights
    if ("weight" %notin% names(popn2)) {
        popn2[, weight := 1]
        popn2[trial == 1, weight := rlnorm(.N, 3.430757, 0.2899917)]
        popn2[trial == 2, weight := rlnorm(.N, 4.457181, 0.3330279)]
    }
    
    popn2
}
