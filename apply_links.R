#' Ensure linked traits are correctly applied
#'
#' @description Take `params` and ensure that linked traits are correct by
#'   setting `Sigma_G`, `priors`, and `fe_vals.`
#'
#' @param params The list of model parameters
#' @returns A new modified model parameters

source("map_idxs.R")

apply_links <- function(params) {
    {
        model_traits <- params$model_traits
        IEs <- params$IEs
        FEs <- params$FEs
        weight_is_nested <- params$weight_is_nested
        fe_vals <- params$fe_vals
        Sigma_G <- params$Sigma_G
        Sigma_E <- params$Sigma_E
        cov_G <- params$cov_G
        cov_E <- params$cov_E
    }


    # Make explicit copies because data.tables work by reference
    params2 <- copy(params)

    # These should start at 0
    params2$fe_vals <- fe_vals * 0
    params2$Sigma_G <- Sigma_G * 0
    params2$Sigma_E <- Sigma_E * 0
    params2$cov_G <- cov_G * 0
    params2$cov_E <- cov_E * 0

    # All changes made to priors also happen in params2 (since DT works by
    # reference)
    priors <- params2$priors

    # Get the first letter of each model traitname for convenience
    t1 <- names(model_traits)
    ntraits <- length(model_traits)

    # Turn "sit", "sittt", or "si__t" etc. into a map
    get_idxs <- function(IEs) {
        s <- str_chars(IEs)
        i1 <- match(s, t1) |> discard(is.na)
        i2 <- if (length(s) == ntraits) {
            match(t1[str_detect(s, "_", negate = TRUE)], t1)
        } else {
            match(s, t1)
        }
        list(i1 = i1, i2 = i2)
    }

    # Individual Effects ----
    if (IEs %notin% c("", "none", NA)) {
        i1 <- get_lidx(IEs, model_traits)
        i2 <- get_idx(IEs, model_traits)

        # Copy rows and columns in Sigma and cov matrix
        params2$Sigma_G[i2, i2] <- Sigma_G[i1, i1, drop = FALSE]
        params2$Sigma_E[i2, i2] <- Sigma_E[i1, i1, drop = FALSE]
        params2$cov_G[i2, i2] <- cov_G[i1, i1, drop = FALSE]
        params2$cov_E[i2, i2] <- cov_E[i1, i1, drop = FALSE]

        # Be careful here in case any traits aren't used
        set_diag_to_cov <- function(x) {
            diag(x) <- str_replace(diag(x), "r", "cov")
            x
        }
        get_LT_pars <- function(x) c(diag(x), t(x)[lower.tri(x)])
        duplicate_GE <- function(x) c(x, str_replace(x, "G", "E"))

        pars <- expand.grid("r_G_", t1, t1) |>
            apply(1, str_flatten) |>
            matrix(ntraits, ntraits) |>
            set_diag_to_cov() |>
            get_LT_pars() |>
            duplicate_GE()

        parvals <- data.table(par = pars,
                              val = c(diag(Sigma_G),
                                      Sigma_G[lower.tri(Sigma_G)],
                                      diag(Sigma_E),
                                      Sigma_E[lower.tri(Sigma_E)]))

        priors[parameter %in% parvals$par,
               true_val := fifelse(use, parvals$val, true_val)]
    }

    get_priors <- function(fe = "trial") {
        if (fe %notin% names(FEs)) return("XYZ")
        FEs[[fe]] |> str_chars() |> unique() |> str_subset("_", negate = TRUE) |>
            str_c(fe, "_", x = _)
    }

    # Fixed Effects ----
    walk(names(FEs), \(FE) {
        i1 <- get_lidx(FEs[[FE]], model_traits)
        i2 <- get_idx(FEs[[FE]], model_traits)

        params2$fe_vals[FE, i2] <<- params$fe_vals[FE, i1] |> unname()

        used <- get_priors(FE)
        priors[str_starts(parameter, str_c(FE, "_")),
               use := parameter %in% used]
        priors[parameter %in% used,
               true_val := params$fe_vals[FE, unique(i1)] |> unname()]
    })

    # Ensure weight_is_nested is respected
    wt <- if (weight_is_nested) "weight\\d_" else "weight_"
    priors[str_starts(parameter, wt), use := FALSE]

    # Copy updated values back into params2
    params2$fe_vals <- fe_vals

    params2
}

