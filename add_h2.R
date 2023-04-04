add_h2 <- function(x, pars) {
    
    if (all(c("cov_G_ss", "cov_E_ss") %in% names(x))) {
        x[, h2_ss := 1 / (1 + cov_E_ss / cov_G_ss)]
        pars <- c(pars, "h2_ss")
        
    }
    
    if (all(c("cov_G_ii", "cov_E_ii") %in% names(x))) {
        x[, h2_ii := 1 / (1 + cov_E_ii / cov_G_ii)]
        pars <- c(pars, "h2_ii")
    }
    
    if (all(c("cov_G_rr", "cov_E_rr") %in% names(x))) {
        x[, h2_rr := 1 / (1 + cov_E_rr / cov_G_rr)]
        pars <- c(pars, "h2_rr")
    }
    
    pars
}