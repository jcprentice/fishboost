set_weights <- function(popn, params) {
    if ("weight" %in% names(popn)) {
        return(popn)
    }
    
    # If weight is missing, generate new values by sampling directly from the FB
    # data. I'd like to generate fresh values from a log-normal distribution,
    # but this seems to result in a much large range of values than I want.
    
    popn[, weight := 1]
    # popn[trial == 1, weight := rlnorm(.N, 3.430757, 0.2899917)]
    # popn[trial == 2, weight := rlnorm(.N, 4.457181, 0.3330279)]

    x <- readRDS("fb_data/fb_12_rpw.rds")
    popn[trial == 1, weight := sample(x[trial == 1, weight], .N, replace = TRUE)]
    popn[trial == 2, weight := sample(x[trial == 2, weight], .N, replace = TRUE)]
    setcolorder(popn, "weight", after = "group")
    
    popn
}
