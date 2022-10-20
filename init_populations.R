# Initialise populations
# Does some trickery to choose the right functions to call
init_pop <- function(traits, params) {
    compartments <- params$compartments
    timings      <- params$timings
    traitnames   <- params$traitnames

    # only want to work with the progeny
    pop <- traits[sdp == "progeny"]

    # add extra columns
    pop[, status := factor("S", levels = compartments)]
    pop[, (timings) := NA_real_]
    pop[, group_inf := 0.0]
    pop[, generation := NA_integer_]

    # handle donors
    pop[donor == 1L, c("generation", "infected_by") := list(1L, 0L)]

    # select appropriate path generating function
    gp <- get(paste0("generate_", tolower(params$model_type), "_path"))

    for (idx in pop[, .I[donor == 1L]]) {
        set(pop, idx, c("status", timings), gp(0.0, pop, idx, params))
    }

    pop
}

