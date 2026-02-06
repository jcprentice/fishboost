set_ge_opts <- function(params) {
    msgs    <- params$msgs
    ge_opts <- params$ge_opts
    
    if (is.null(ge_opts) || ge_opts == "" || length(ge_opts) == 0) {
        return(params)
    }
    
    params2 <- copy(params)
    
    if (any(str_starts(ge_opts, "no_ev"))) {
        evs <- ge_opts |> str_subset("no_ev") |> str_remove("no_ev_") |> str_chars()
        
        if (msgs) message("Removing VarE from ", evs)
        
        walk(evs, \(x) {
            sit <- c("s", "i", "t")
            cpattern <- str_glue("cov_E_{x}{sit}|r_E_{x}{sit}|r_E_{sit}{x}") |>
                str_c(collapse = "|")
            params2$priors[str_detect(parameter, cpattern), use := FALSE]
        })
        params2$h2[] <- 0
    }
    
    if ("gt_only" %in% ge_opts) {
        if (msgs) message("Removing all VarE")
        
        params2$priors[str_detect(parameter, "cov_E|r_E"), use := FALSE]
        params2$h2[] <- 0
    }
    
    if ("pt_only" %in% ge_opts) {
        if (msgs) message("Removing all VarG")
        
        params2$priors[str_detect(parameter, "cov_G|r_G"), use := FALSE]
        params2$h2[] <- 0
    }
    
    if ("e1" %in% ge_opts) {
        if (msgs) message("Setting all environmental covariance to I")
            
        params2$priors[str_detect(parameter, "cov_E"), `:=`(type = "Fixed", val1 = 1, val2 = 1)]
        params2$priors[str_detect(parameter, "r_E"),   `:=`(type = "Fixed", val1 = 0, val2 = 0)]
        params2$cor_E[] <- with(params2, diag(1, nrow(cor_E), ncol(cor_E)))
        params2$h2 <- with(params2, Sigma_G / (Sigma_G + Sigma_E))
    }
    
    params2
}
