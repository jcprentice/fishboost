# Initialise populations
# Does some trickery to choose the right functions to call
init_popn <- function(popn, params) {
    compartments <- params$compartments
    timings      <- params$timings
    traitnames   <- params$traitnames
    model_type   <- params$model_type
    
    # Only want to work with the progeny
    popn2 <- popn[sdp == "progeny"]
    
    # Add extra columns
    popn2[, (timings) := NA_real_]
    popn2[, `:=`(status = factor("S", levels = compartments),
                 group_inf = 0.0,
                 generation = NA_integer_)]
    
    # Handle donors
    popn2[donor == 1L, `:=`(generation = 1L,
                            infected_by = 0L)]
    
    # Select appropriate path generating function
    gp <- get(str_c("generate_", str_to_lower(model_type), "_path"))
    
    # Generate paths for donors
    walk(popn2[, .I[donor == 1L]], \(i) {
        set(popn2, i, c("status", timings), gp(0.0, popn2, i, params))
    })
    
    popn2
}

