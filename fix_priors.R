fix_priors(params) {
    pars <- params$priors$parameter
    
    fe_pars  <- pars[grepl("^trial|^donor|^txd", pars)]
    cov_pars <- pars[grepl("^cov", pars)]
    r_pars   <- pars[grepl("r_G|r_E", pars)]
    
    other_pars <- setdiff(pars, c(fe_pars, cov_pars, r_pars))
    
    # Ensure true_val is in (val1, val2)
    params$priors[parameter %in% fe_pars,
           `:=`(val1 = pmin(val1, true_val - 2),
                val2 = pmax(val2, true_val - 2))]
    
    # Ensure true_val is in (0, val2)
    params$priors[parameter %in% cov_pars,
                  val2 := pmax(val2, true_val + 1)]
    
    message("Do something here?")
    params$priors[parameter %in% other_pars,
                  val2 := val2] # FIXME
    
    params
}