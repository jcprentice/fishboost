library(stringr)

#' Take BICI parameters and make them consistent with internal parameter names,
#' e.g.
#' - "G_1" -> "Group effect 1"
#' - "\\mu^weight1s" -> "weight1_s"
#' - "\\Omega^gen_sg,sg" -> "cov_G_ss"
#' - "\\omega^gen_sg,ig" -> "r_G_si"
#'
#' @param pars A character vector of BICI parameter names
#'
#' @returns A new character vector of modified parameter names

rename_bici_pars <- function(pars) {
    pars |> str_replace_all(
        c("\\\\" = "",
          "\\^" = "_",
          "cv_G" = "sigma",
          "^G_" = "Group effect ",
          # "infrat" = "Inf ratio",
          # "Ω" = "Omega", "ω" = "omega", "μ" = "mu",
          "mu_weight([12]?)([slidt])" = "weight\\1_\\2",
          "Omega_gen_(.)g,(.)g" = "cov_G_\\1\\2",
          "Omega_env_(.)e,(.)e" = "cov_E_\\1\\2",
          "omega_gen_(.)g,(.)g" = "r_G_\\1\\2",
          "omega_env_(.)e,(.)e" = "r_E_\\1\\2"))
}

#' Rename parameters
#'
#' Turns the short name parameters into something suitable for printing, e.g.
#' - "cov_G_ss" -> "Var G (sus)"
#' - "latent_period" -> "Latent Period"
#' - "trial_i" -> "Trial Infectivity"
#'
#' @param pars A character vector of parameter names
#'
#' @returns A new character vector of modified parameter names

html_names <- function(pars) {
    pars |> str_replace_all(
        c("latent_period" = "LP",
          "detection_period" = "DP",
          "removal_period" = "RP",
          "Tr(.),Don" = "(Tr\\1, Seeder)",
          "Tr(.),Rec" = "(Tr\\1, Contact)",
          "beta_Tr(.)" = "beta (Tr\\1)",
          "infrat" = "Seeder inf ratio",
          "^trial" = "Trial",
          "^donor" = "Donor",
          "^txd" = "Trial x Donor",
          "weight_" = "Weight_",
          "weight(.)_" = "Weight (Tr\\1,_",
          "Tr([12])" = "Trial \\1",
          "sigma" = "Group Effect",
          "^cov_" = "Var_",
          "^r_" = "Cor_",
          "_G_" = "<sub>A</sub>_",
          "_E_" = "<sub>E</sub>_",
          "_P_" = "<sub>P</sub>_",
          "^h2_" = "h<sup>2</sup>_",
          "_s+$" = "_(Sus)",
          "_i+$" = "_(Inf)",
          "_l+$" = "_(Lat)",
          "_d+$" = "_(Det)",
          "_t+$" = "_(End)",
          "_si$" = "_(Sus, Inf)",
          "_sl$" = "_(Sus, Lat)",
          "_sd$" = "_(Sus, Det)",
          "_st$" = "_(Sus, End)",
          "_il$" = "_(Inf, Lat)",
          "_id$" = "_(Inf, Det)",
          "_it$" = "_(Inf, End)",
          "_ld$" = "_(Lat, Det)",
          "_lt$" = "_(Lat, End)",
          "_dt$" = "_(Det, End)",
          ",_\\(" = ",_",
          "_" = " "))
}

pretty_names <- function(pars) {
    pars |> str_replace_all(
        c("latent_period" = "LP",
          "detection_period" = "DP",
          "removal_period" = "RP",
          "Tr(.),Don" = "(Tr\\1, Seeder)",
          "Tr(.),Rec" = "(Tr\\1, Contact)",
          "beta_Tr(.)" = "beta (Tr\\1)",
          "infrat" = "Seeder inf ratio",
          "^trial" = "Trial",
          "^donor" = "Donor",
          "^txd" = "Trial x Donor",
          "weight_" = "Weight_",
          "weight(.)_" = "Weight (Tr\\1,_",
          "Tr([12])" = "Trial \\1",
          "sigma" = "Group Effect",
          "^cov_" = "Var_",
          "^r_" = "Cor_",
          "_G_" = "_Gen_",
          "_E_" = "_Env_",
          "_P_" = "_PT_",
          "_s+$" = "_(Sus)",
          "_i+$" = "_(Inf)",
          "_l+$" = "_(Lat)",
          "_d+$" = "_(Det)",
          "_t+$" = "_(End)",
          "_si$" = "_(Sus, Inf)",
          "_sl$" = "_(Sus, Lat)",
          "_sd$" = "_(Sus, Det)",
          "_st$" = "_(Sus, End)",
          "_il$" = "_(Inf, Lat)",
          "_id$" = "_(Inf, Det)",
          "_it$" = "_(Inf, End)",
          "_ld$" = "_(Lat, Det)",
          "_lt$" = "_(Lat, End)",
          "_dt$" = "_(Det, End)",
          ",_\\(" = ",_",
          "_" = " "))
}

sildt <- c("Susceptibility", "Infectivity", "Latency", "Detectability", "Tolerance")
sildt1 <- c("s", "i", "l", "d", "t")
sildt2 <- c("ss", "ii", "ll", "dd", "tt")
xy <- c("si", "sl", "sd", "st", "il", "id", "it", "ld", "lt", "dt")

fes <- expand.grid(c("s", "i", "l", "d", "t"),
                   c("trial", "donor", "txd", "weight", "weight1", "weight2")) |>
    rev() |> apply(1, str_flatten, "_")

# Preferred parameter order

param_order <- c(
    str_c("cov_G_", sildt2), str_c("r_G_", xy),
    str_c("cov_E_", sildt2), str_c("r_E_", xy),
    str_c("cov_P_", sildt2), str_c("h2_", sildt2),
    "beta_Tr1", "beta_Tr2", "sigma", "infrat",
    "LP_Tr1,Don", "LP_Tr1,Rec", "LP_Tr2,Don", "LP_Tr2,Rec",
    "DP_Tr1,Don", "DP_Tr1,Rec", "DP_Tr2,Don", "DP_Tr2,Rec",
    "RP_Tr1,Don", "RP_Tr1,Rec", "RP_Tr2,Don", "RP_Tr2,Rec",
    fes, "MVPSF"
) |>
    str_replace_all(c("LP" = "latent_period",
                      "DP" = "detection_period",
                      "RP" = "removal_period"))

full_param_order <- pretty_names(param_order)

