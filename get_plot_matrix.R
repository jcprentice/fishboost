{
    library(stringr)
    library(purrr)
    source("utils.R")
}

get_plot_matrix <- function(pars, plt_shape = "traits") {
    pars <- setdiff(pars, "title_plt")

    trials <- str_extract(pars, "beta_Tr(.)", group = 1) |>
        discard(is.na) |> unique() |> as.integer() |> sort()

    plt_mat <- if (plt_shape == "traits") {

        sildt1 <- str_chars("sildt")
        sildt2 <- str_c(sildt1, sildt1)
        any_non_empty <- function(x) any(x != "empty")

        cov_pars <- c(str_c("cov_G_", sildt2),
                      "r_G_si", "r_G_st", "empty", "empty", "r_G_it",
                      str_c("cov_E_", sildt2))
        # cov_pars <- c(cov_pars, str_c("cov_P_", sildt2), str_c("h2_", sildt2))

        model_pars <- c(
            "sigma",  "beta_Tr1", "LP_Tr1,Don", "DP_Tr1,Don", "RP_Tr1,Don",
            "infrat", "empty",    "LP_Tr1,Rec", "DP_Tr1,Rec", "RP_Tr1,Rec",
            "sigma",  "beta_Tr2", "LP_Tr2,Don", "DP_Tr2,Don", "RP_Tr2,Don",
            "infrat", "empty",    "LP_Tr2,Rec", "DP_Tr2,Rec", "RP_Tr2,Rec")

        # Remove repeated sigma and infrat if only one trial
        if (identical(trials, 1L)) {
            model_pars[c(11, 16)] <- "empty"
        } else if (identical(trials, 2L)) {
            model_pars[c(1, 6)] <- "empty"
        }

        fes <- expand.grid(sildt1,
                           c("trial", "donor", "txd", "weight", "weight1", "weight2")) |>
            rev() |> apply(1, str_flatten, "_")

        plt_names <- c(cov_pars, model_pars, fes)

        # Some entries like "trial_s" might be missing
        plt_names[plt_names %notin% pars] <- "empty"

        # This clips any rows or columns that are entirely empty
        plt_mat <- matrix(plt_names, ncol = 5, byrow = TRUE)

        plt_mat[which(apply(plt_mat, 1, any_non_empty)),
                which(apply(plt_mat, 2, any_non_empty))]

    } else if (plt_shape == "compact") {
        plt_names <- c(
            "cov_G_ss", "cov_G_ii", "cov_G_tt", "r_G_si", "r_G_st", "r_G_it",
            "cov_E_ss", "cov_E_ii", "cov_E_tt", "r_E_si", "r_E_st", "r_E_it",
            if (identical(trials, 1L)) {
                c("LP_Tr1,Don", "DP_Tr1,Don", "RP_Tr1,Don", "weight1_s", "weight1_i", "weight1_t",
                  "LP_Tr1,Rec", "DP_Tr1,Rec", "RP_Tr1,Rec", "beta_Tr1",  "infrat",    "sigma")
            } else if (identical(trials, 2L)) {
                c("LP_Tr2,Don", "DP_Tr2,Don", "RP_Tr2,Don", "weight2_s", "weight2_i", "weight2_t",
                  "LP_Tr2,Rec", "DP_Tr2,Rec", "RP_Tr2,Rec", "beta_Tr2",  "infrat",    "sigma")
            } else {
                c("LP_Tr1,Don", "DP_Tr1,Don", "RP_Tr1,Don", "weight1_s", "weight1_i", "weight1_t",
                  "LP_Tr1,Rec", "DP_Tr1,Rec", "RP_Tr1,Rec", "weight2_s", "weight2_i", "weight2_t",
                  "LP_Tr2,Don", "DP_Tr2,Don", "RP_Tr2,Don", "beta_Tr1",  "beta_Tr2",  "infrat",
                  "LP_Tr2,Rec", "DP_Tr2,Rec", "RP_Tr2,Rec", "sigma",     "empty",     "empty")
            }
        )

        plt_names[plt_names %notin% pars] <- "empty"

        matrix(plt_names, ncol = 6, byrow = TRUE)

    } else {
        # plt_shape == "rectangle"
        nr <- ceiling(sqrt(length(pars)))
        nc <- ceiling(length(pars) / nr)

        matrix(c(pars, rep("empty", nr * nc - length(pars))),
               nrow = nr,
               ncol = nc)
    }

    list(plt_names = c(t(plt_mat)),
         nc = ncol(plt_mat),
         nr = nrow(plt_mat))
}
