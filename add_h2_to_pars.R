library(data.table)
library(stringr)

#' Add heritability and phenotypes
#'
#' @description Modify a data.table `x` of posterior values to add new phenotypic
#' variance and heritability columns.
#'
#' @param x A data.table of posterior values
#' @param pars A character vector of parameters
#'
#' @returns A new pars vector

add_h2_to_pars <- function(x, pars) {

    y <- data.table(expand.grid(sit = c("s", "i", "t"), ae = c("a", "e")))
    y[, `:=`(bici = str_c("\\Omega^", sit, ae, ",", sit, ae),
             sire = str_c("cov_", fifelse(ae == "a", "G", "E"), "_", sit, sit))]

    # Prefer working with SIRE's namning scheme
    
    setnames(x, y$bici, y$sire, skip_absent = TRUE)

    if (all(c("cov_G_ss", "cov_E_ss") %in% names(x))) {
        x[, cov_P_ss := cov_G_ss + cov_E_ss]
        x[, h2_ss := cov_G_ss / cov_P_ss]
        pars <- c(pars, "cov_P_ss", "h2_ss")
    }

    if (all(c("cov_G_ii", "cov_E_ii") %in% names(x))) {
        x[, cov_P_ii := cov_G_ii + cov_E_ii]
        x[, h2_ii := cov_G_ii / cov_P_ii]
        pars <- c(pars, "cov_P_ii", "h2_ii")
    }

    if (all(c("cov_G_tt", "cov_E_tt") %in% names(x))) {
        x[, cov_P_tt := cov_G_tt + cov_E_tt]
        x[, h2_tt := cov_G_tt / cov_P_tt]
        pars <- c(pars, "cov_P_tt", "h2_tt")
    }

    pars
}

