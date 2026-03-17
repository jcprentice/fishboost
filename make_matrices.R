# Construct cov_mat and Sigma matrices using var and cor components
make_matrices <- function(model_traits = c("sus", "inf", "lat", "det", "tol"),
                          vars = 1,
                          cors = 0.2) {
    if (FALSE) {
        model_traits <- c("sus", "inf", "lat", "det", "tol")
        vars <- params$vars; cors <- params$cors
        vars <- 1; cors <- 0.2
        vars <- list(sus = 1, inf = 1.5, tol = 0.5, default = 0.2)
        cors <- list(si = -0.3, st = -0.3, it =  0.3, default = 0.15)
    }

    traits1 <- names(model_traits)
    N <- length(model_traits)

    vars <- unlist(vars)
    cors <- unlist(cors)

    Sigma <- matrix(0, N, N, dimnames = list(model_traits, model_traits))

    # First set all values to whatever the default cors is
    if (is.null(names(cors))) {
        Sigma[] <- cors[[1]]
    } else if ("default" %in% names(cors)) {
        Sigma[] <- cors[["default"]]
    }

    # Next set the default vars along the diagonal
    if (is.null(names(vars))) {
        diag(Sigma) <- vars[[1]]
    } else if ("default" %in% names(vars)) {
        diag(Sigma) <- vars[["default"]]
    }

    # Now fill in named values for vars
    for (i in model_traits) {
        if (i %in% names(vars)) {
            Sigma[i, i] <- vars[[i]]
        }
    }

    # Finally fill in named values for cors
    mat <- expand.grid(traits1, traits1) |> rev() |> apply(1, str_flatten) |> matrix(N, N)
    ij <- matrix(1:N^2, N, N)[lower.tri(mat)]
    cor_names <- mat[ij]

    for (i in ij) {
        if (mat[i] %in% names(cors)) {
            Sigma[i] <- cors[[mat[i]]]
        }
    }

    # Convert Sigma into a covariance matrix
    D <- 0 * Sigma
    diag(D) <- sqrt(diag(Sigma))
    R <- Sigma
    diag(R) <- 1
    cov_mat <- D %*% R %*% D

    list(Sigma = Sigma, cov = cov_mat, cor_names = cor_names)
}


# Construct cov_mat and Sigma matrices from the priors DT
make_matrices_from_priors <- function(priors) {
    all_traits <- c("sus", "inf", "lat", "det", "tol")

    xx <- c("ss", "ii", "ll", "dd", "tt")
    xy <- c("si", "sl", "sd", "st", "il", "id", "it", "ld", "lt", "dt")
    covs <- c(str_c("cov_G_", xx), str_c("r_G_", xy),
              str_c("cov_E_", xx), str_c("r_E_", xy))
    walk(covs, \(i) priors[[i]] <<- priors[[i]] %||% 0)

    Sigma_G <- with(priors,
                    matrix(c(cov_G_ss, r_G_si,   r_G_sl,   r_G_sd,   r_G_st,
                             r_G_si,   cov_G_ii, r_G_il,   r_G_id,   r_G_it,
                             r_G_sl,   r_G_il,   cov_G_ll, r_G_ld,   r_G_lt,
                             r_G_sd,   r_G_id,   r_G_ld,   cov_G_dd, r_G_dt,
                             r_G_st,   r_G_it,   r_G_lt,   r_G_dt,   cov_G_tt),
                           5, 5,
                           dimnames = list(all_traits, all_traits)))

    Sigma_E <- with(priors,
                    matrix(c(cov_G_ss, r_G_si,   r_G_sl,   r_G_sd,   r_G_st,
                             r_G_si,   cov_G_ii, r_G_il,   r_G_id,   r_G_it,
                             r_G_sl,   r_G_il,   cov_G_ll, r_G_ld,   r_G_lt,
                             r_G_sd,   r_G_id,   r_G_ld,   cov_G_dd, r_G_dt,
                             r_G_st,   r_G_it,   r_G_lt,   r_G_dt,   cov_G_tt),
                           5, 5,
                           dimnames = list(all_traits, all_traits)))

    D <- 0 * Sigma_G
    diag(D) <- sqrt(diag(Sigma_G))
    R <- Sigma_G
    diag(R) <- 1
    cov_G <- D %*% R %*% D

    D <- 0 * Sigma_E
    diag(D) <- sqrt(diag(Sigma_E))
    R <- Sigma_E
    diag(R) <- 1
    cov_E <- D %*% R %*% D

    list(Sigma_G = Sigma_G,
         Sigma_E = Sigma_E,
         cov_G = cov_G,
         cov_E = cov_E)
}
