# Turn a Sigma matrix into a cov matrix
sigma2cov <- function(Sigma) {
    D <- 0 * Sigma
    diag(D) <- sqrt(diag(Sigma))
    R <- Sigma
    diag(R) <- 1

    D %*% R %*% D
}

# Construct cov_mat and Sigma matrices using var and cor components
make_matrices <- function(model_traits,
                          use_traits = "sit",
                          h2 = 0.5,
                          vars = 1,
                          cors = 0.2) {
    if (FALSE) {
        model_traits <- c(s = "sus", i = "inf", l = "lat", d = "det", t = "tol")
        use_traits <- "sit"
        h2 <- 0.5
        h2 <- params$h2; vars <- params$vars; cors <- params$cors
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

    Sigma_G <- Sigma
    Sigma_E <- Sigma
    diag(Sigma_G) <- diag(Sigma) * 2 * h2
    diag(Sigma_E) <- diag(Sigma) * 2 * (1 - h2)


    # Finally fill in named values for cors
    mat <- expand.grid(traits1, traits1) |> rev() |> apply(1, str_flatten) |> matrix(N, N)
    ij <- matrix(1:N^2, N, N)[lower.tri(mat)]
    cor_names <- mat[ij]

    for (i in ij) {
        if (mat[[i]] %in% names(cors)) {
            Sigma[[i]] <- cors[[mat[[i]]]]
        }
    }

    # Apply use_traits (zero out unused traits)
    unused_traits <- model_traits[setdiff(traits1, str_chars(use_traits))]
    for (i in unused_traits) {
        Sigma[i, ] <- 0
        Sigma[, i] <- 0
    }

    list(Sigma_G = Sigma_G,
         Sigma_E = Sigma_E,
         cov_G = sigma2cov(Sigma_G),
         cov_E = sigma2cov(Sigma_E),
         cor_names = cor_names)
}


# Construct cov_mat and Sigma matrices from the priors DT
make_matrices_from_priors <- function(priors1, params) {
    all_traits <- c("sus", "inf", "lat", "det", "tol")
    N <- length(all_traits)

    xx <- c("ss", "ii", "ll", "dd", "tt")
    xy <- c("si", "sl", "sd", "st", "il", "id", "it", "ld", "lt", "dt")
    covs <- c(str_c("cov_G_", xx), str_c("r_G_", xy),
              str_c("cov_E_", xx), str_c("r_E_", xy))
    walk(covs, \(i) priors1[[i]] <<- priors1[[i]] %||% 0)

    Sigma_G <- with(priors1, matrix(
        c(cov_G_ss, r_G_si,   r_G_sl,   r_G_sd,   r_G_st,
          r_G_si,   cov_G_ii, r_G_il,   r_G_id,   r_G_it,
          r_G_sl,   r_G_il,   cov_G_ll, r_G_ld,   r_G_lt,
          r_G_sd,   r_G_id,   r_G_ld,   cov_G_dd, r_G_dt,
          r_G_st,   r_G_it,   r_G_lt,   r_G_dt,   cov_G_tt),
        N, N,
        dimnames = list(all_traits, all_traits)
    ))

    Sigma_E <- with(priors1, matrix(
        c(cov_E_ss, r_E_si,   r_E_sl,   r_E_sd,   r_E_st,
          r_E_si,   cov_E_ii, r_E_il,   r_E_id,   r_E_it,
          r_E_sl,   r_E_il,   cov_E_ll, r_E_ld,   r_E_lt,
          r_E_sd,   r_E_id,   r_E_ld,   cov_E_dd, r_E_dt,
          r_E_st,   r_E_it,   r_E_lt,   r_E_dt,   cov_E_tt),
        N, N,
        dimnames = list(all_traits, all_traits)
    ))

    if (params$override_h2 %||% FALSE) {
        h2 <- params$h2
        message("- Overriding with h2 = ", h2)
        vars <- diag(Sigma_G) + diag(Sigma_E)
        diag(Sigma_G) <- diag(Sigma_G) * vars * h2
        diag(Sigma_E) <- diag(Sigma_E) * vars * (1 - h2)

        if (sum(diag(Sigma_G)) < 1e-3) {
            diag(Sigma_G) <- 1e-3
        }
        if (sum(diag(Sigma_E)) < 1e-3) {
            diag(Sigma_E) <- 1e-3
        }
    }

    list(Sigma_G = Sigma_G,
         Sigma_E = Sigma_E,
         cov_G = sigma2cov(Sigma_G),
         cov_E = sigma2cov(Sigma_E))
}
