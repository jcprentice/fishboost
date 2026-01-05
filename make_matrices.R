# Construct cov_mat and Sigma matrices using var and cor components
make_matrices <- function(model_traits = c("sus", "inf", "lat", "det", "tol"),
                          use_traits = "sildt",
                          vars = 1, cors = 0.2) {
    if (FALSE) {
        model_traits <- params$model_traits
        vars <- params$vars
        cors <- params$cors
    }
    all_traits <-  c("sus", "inf", "lat", "det", "tol")
    sildt <- c("s", "i", "l", "d", "t")
    idxs <- match(model_traits, all_traits)
    
    vars <- unlist(vars)
    cors <- unlist(cors)
    
    if (FALSE) {
        vars <- c(sus = 1, inf = 1.5, lat = 0, det = 0, tol = 0.5)
        cors <- c(si = -0.3, sl = 0.2, sd = 0.2, st = -0.3,
                  il = 0.2, id = 0.2, it =  0.3,
                  ld = 0.2, lt =  0.2,
                  dt =  0.2)
    }
    if (FALSE) {
        vars <- 1; cors <- 0.2
    }
    
    # Get the matrix idxs to assign cors to
    mat <- matrix(1:25, 5, 5, dimnames = list(all_traits, all_traits))
    mat_idxs <- mat[idxs, idxs][lower.tri(mat[idxs, idxs])]
    cns <- expand.grid(sildt, sildt) |> rev()  |> apply(1, str_flatten) |> _[mat_idxs]
    
    # Create Sigma and the correlation matrix with `cors`
    D <- cor_mat <- matrix(0, 5L, 5L, dimnames = list(all_traits, all_traits))
    if (is.null(names(cors))) {
        if (length(cors) %notin% c(1, length(mat_idxs))) {
            capture_message(cors)
            stop(str_glue("length(cors) = {lc}, expecting 1, {lmi}, or named vector",
                          lc = length(cors), lmi = length(mat_idxs)))
        }
        cor_mat[mat_idxs] <- cors
    } else {
        expand.grid(sildt, sildt) |>
            rev() |> _[mat_idxs, ] |>
            pwalk(~ {
                cn <- c(.x, .y) |> str_flatten()
                if (cn %in% names(cors)) {
                    cor_mat[.x, .y] <<- cors[[cn]]
                }
            })
    }    
    cor_mat <- cor_mat + t(cor_mat)
    diag(cor_mat) <- 1
    
    # Sigma starts off as cor_mat, but overwrite with the variances
    Sigma <- cor_mat
    diag(Sigma) <- 0
    if (is.null(names(vars))) {
        if (length(vars) %notin% c(1, length(idxs))) {
            capture_message(vars)
            stop(str_glue("length(vars) = {lv}, expecting 1, {li}, or named vector",
                          lv = length(vars), li = length(idxs)))
        }
        diag(Sigma)[idxs] <- vars
    } else {
        walk(names(vars), ~ {Sigma[.x, .x] <<- vars[[.x]]})
    }
    
    # D is the std dev, or the sqrt of the diag of Sigma
    diag(D) <- sqrt(diag(Sigma))
    
    # Use `vars` to make the covariance matrix
    cov_mat <- D %*% cor_mat %*% D
    
    # Remove unused traits
    idxs0 <- str_which(use_traits, str_1st(all_traits), negate = TRUE) |>
        union(which(diag(Sigma) == 0))
    Sigma[idxs0, ] <- 0
    Sigma[, idxs0] <- 0
    cov_mat[idxs0, ] <- 0
    cov_mat[, idxs0] <- 0
    
    # Discard rows and columns not in model_traits
    Sigma <- Sigma[idxs, idxs]
    cov_mat <- cov_mat[idxs, idxs]
    
    list(Sigma = Sigma, cov = cov_mat, cns = cns)
}


# Construct cov_mat and Sigma matrices from the priors DT
make_matrices_from_priors <- function(priors) {
    all_traits <- c("sus", "inf", "lat", "det", "tol")
    
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
    
    D <- diag(sqrt(diag(Sigma_G)))
    R <- Sigma_G
    diag(R) <- 1
    cov_G <- D %*% R %*% D
    dimnames(cov_G) <- list(all_traits, all_traits)
    
    D <- diag(sqrt(diag(Sigma_E)))
    R <- Sigma_E
    diag(R) <- 1
    cov_E <- D %*% R %*% D
    dimnames(cov_E) <- list(all_traits, all_traits)
    
    list(Sigma_G = Sigma_G,
         Sigma_E = Sigma_E,
         cov_G = cov_G,
         cov_E = cov_E)
}
