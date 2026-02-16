library(stringr)

# Turn BICI pars into something useful
rename_bici_pars <- function(pars) {
    pars |> str_replace_all(
        c("\\\\" = "",
          "\\^" = "_",
          "cv_G" = "sigma",
          "^G_" = "Group effect ",
          # "infrat" = "Inf ratio",
          "mu_weight([12]?)([slidt])" = "weight\\1_\\2",
          "Omega_gen_(.)g,(.)g" = "cov_G_\\1\\2",
          "Omega_env_(.)e,(.)e" = "cov_E_\\1\\2",
          "omega_gen_(.)g,(.)g" = "r_G_\\1\\2",
          "omega_env_(.)e,(.)e" = "r_E_\\1\\2"))
}

# This turns the short name into something suitable for printing
rename_pars <- function(pars) {
    pars |> str_replace_all(
        c("latent_period" = "Latent Period (days)",
          "LP_shape" = "LP shape",
          "detection_period" = "Detection Period (days)",
          "DP_shape" = "DP shape",
          "removal_period" = "Removal Period (days)",
          "RP_shape" = "RP shape",
          "infrat" = "Inf ratio",
          "_s$" = "_Susceptibility",
          "_l$" = "_Latency",
          "_i$" = "_Infectivity",
          "_d$" = "_Detectability",
          "_t$" = "_Tolerance",
          "^trial" = "Trial",
          "^donor" = "Donor",
          "^txd" = "TxD",
          "sigma" = "Group Effect",
          "^cov_" = "Var_",
          "^r_" = "Cor_",
          "^h2_" = "h\\^2_",
          "_ss$" = "_(Sus)",
          "_ii$" = "_(Inf)",
          "_tt$" = "_(Tol)",
          "_si$" = "_(Sus, Inf)",
          "_st$" = "_(Sus, Tol)",
          "_it$" = "_(Inf, Tol)",
          "_" = " "))
}

sildt <- c("Susceptibility", "Infectivity", "Latency", "Detectability", "Tolerance")
sildt1 <- sildt |> str_sub(1, 1) |> str_to_lower()

fes <- expand.grid(sildt1,
                   c("trial", "donor", "txd", "weight", "weight1", "weight2")) |>
    rev() |> apply(1, str_flatten, "_")

# Preferred parameter order
param_order <- c(
    # BICI parameters
    "beta_Tr1", "beta_Tr2", "LP_Tr1,Don", "LP_Tr1,Rec", "LP_Tr2,Don",
    "LP_Tr2,Rec", "DP_Tr1,Don", "DP_Tr1,Rec", "DP_Tr2,Don", "DP_Tr2,Rec",
    "RP_Tr1,Don", "RP_Tr1,Rec", "RP_Tr2,Don", "RP_Tr2,Rec", "gamma_shape_Tr1,Don",
    "gamma_shape_Tr1,Rec", "gamma_shape_Tr2,Don", "gamma_shape_Tr2,Rec",

    # SIRE parameters
    "beta", "latent_period", "LP_shape", "detection_period", "DP_shape",
    "removal_period", "RP_shape", "sigma",
    "cov_G_ss", "cov_G_ii",  "cov_G_tt", "r_G_si", "r_G_st", "r_G_it",
    "cov_E_ss", "cov_E_ii",  "cov_E_tt", "r_E_si", "r_E_st", "r_E_it",
    "cov_P_ss", "cov_P_ii",  "cov_P_tt", "h2_ss", "h2_ii", "h2_tt",
    fes,
    "MVPSF")

full_param_order <- rename_pars(param_order)

