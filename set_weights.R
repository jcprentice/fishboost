set_weights <- function(popn, params) {
    if ("weight" %in% names(popn)) {
        return(popn)
    }
    
    # If weights aren't assigned, assign new ones
    
    weight_type <- params$weight_type %||% "none"
    
    popn[, weight := 1]
    setcolorder(popn, "weight", after = "group")
    
    if (weight_type == "sampled") {
        fb <- readRDS("fb_data/fb_12_rpw.rds")
        popn[trial == 1, weight := sample(fb[trial == 1, weight], .N, replace = TRUE)]
        popn[trial == 2, weight := sample(fb[trial == 2, weight], .N, replace = TRUE)]
    } else {
        # Generate a truncated log-normal distribution
        fb <- readRDS("fb_data/fb_12_rpw.rds")
        
        trials <- popn[sdp == "progeny", unique(trial)]
        
        walk(trials, \(tr) {
            fit <- fb[trial == tr, MASS::fitdistr(weight, "log-normal")]
            rng <- fb[trial == tr, range(weight)]
            wts <- seq(rng[[1]], rng[[2]], 0.1)
            popn[trial == tr,
                 weight := sample(wts, .N, replace = TRUE,
                                  prob = dlnorm(wts,
                                                fit$estimate[[1]],
                                                fit$estimate[[2]]))]
        })
    }
    
    popn
}
