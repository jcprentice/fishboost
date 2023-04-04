# Find the no. of individuals capable of infecting susceptibles at some point
# (so includes those currently exposed and not yet infectious).
get_seir_res_infectives <- function(pop) {
    pop[, sum(status %in% c("E", "I"))]
}


# A Susceptible individual's future disease trajectory is fixed at the point of
# exposure.
generate_seir_res_path <- function(epi_time, pop, id, params) {
    with(params, {
        Tinf <- epi_time
        Tsym <- Tinf + rgamma(1L, r_eta_shape, r_eta_rate) * pop$latency[id]
        Trec <- Tsym + rgamma(1L, r_gamma_shape, r_gamma_rate) * pop$recoverability[id]

        # return
        list("E", Tinf, Tsym, Trec)
    })
}


# Find the time of the next non-infection event, and the id of the individual
next_seir_res_ni_event <- function(pop, epi_time) {
    as.list(pop[, .(.I,
                    Tsym2 = fifelse(Tsym > epi_time, Tsym, Inf),
                    Trec2 = fifelse(Trec > epi_time, Trec, Inf))
    ][, .(I, Tmin = pmin(Tsym2, Trec2))
    ][, list(t_next_event = min(Tmin, na.rm = TRUE), id_next_event = which.min(Tmin))])
}


# Main model ----
model_SEIR_res <- function(traits, params) {
    message("Simulating an SEIR_res epidemic ...")

    # copy necessary parameters
    r_beta   <- params$r_beta
    r_lambda <- 0.1 #params$r_lambda
    DEBUG    <- params$DEBUG

    # initialise populations ----
    pop <- init_pop(traits, params)
    # group_inf is now the reservoir
    pop[, group_inf := 0.0]

    # Start epidemic simulation loop
    dt <- 0.0
    epi_time <- 0.0
    while (get_seir_res_infectives(pop) > 0L) {
    # {
        if (DEBUG) message("time = ", signif(epi_time, 5))

        # Calculate infection rates in each group
        # this is the sum of the log infectivitities
        pop[, group_inf := group_inf * exp(- r_lambda * dt) +
                100 * GE * mean(fifelse(status == "I", infectivity, 0.0)) * r_beta * dt,
            by = group][]
        if (DEBUG) print(pop[, mean(group_inf), by = group])

        # if S, infection at rate beta SI
        pop[, event_rate := fifelse(status == "S", susceptibility * group_inf, 0.0)]

        # id and time of next non-infection event
        ni_event      <- next_seir_res_ni_event(pop, epi_time)
        t_next_event  <- ni_event$t_next_event
        id_next_event <- ni_event$id_next_event

        if (DEBUG) message("next NI event id = ", id_next_event, " at t = ", t_next_event)

        # generate random timestep ----
        total_event_rate <- sum(pop$event_rate)

        # calculate dt if infections event rate > 0
        if (total_event_rate > 0.0) {
            dt <- rexp(1L, rate = total_event_rate)
        } else {
            dt <- 1e6
        }

        if (DEBUG) message("Total Infections Event Rate = ", signif(total_event_rate, 5))

        # check if next event is infection or non-infection
        if (epi_time + dt < t_next_event) {
            if (DEBUG) message("next event is infection at t = ", signif(epi_time, 5))

            epi_time <- epi_time + dt

            # randomly select individual
            id_next_event <- sample(nrow(pop),
                                    size = 1L,
                                    prob = pop$event_rate)

            group_id <- pop$group[id_next_event]
            infectives <- pop[, .(.I, group, status, infectivity)][group == group_id & status == "I"]

            set(pop, id_next_event, c("status", "Tinf", "Tsym", "Trec"),
                generate_seir_res_path(epi_time, pop, id_next_event, params))

            if (DEBUG) message("ID ", id_next_event, ": S -> E")
        } else {
            if (DEBUG) message("next event is non-infection at t = ", signif(t_next_event, 5))

            dt <- t_next_event - epi_time
            epi_time <- t_next_event

            status <- pop[id_next_event, status]

            if (status == "E") {
                if (DEBUG) message("ID ", id_next_event, ": E -> I")
                set(pop, id_next_event, "status", "I")
            } else if (status == "I") {
                if (DEBUG) message("ID ", id_next_event, ": I -> R")
                set(pop, id_next_event, "status", "R")
            } else {
                print(pop[, .(group, donor, status, Tinf, Tsym, Trec, group_inf, event_rate)])
                print(pop[id_next_event, .(group, donor, status, Tinf, Tsym, Trec, group_inf, event_rate)])
                stop(paste0("selected ID ", id_next_event, "... unexpected event!"))
                break
            }
        }

        if (DEBUG) {
            print(pop[, .(group, donor, status, Tinf, Tsym, Trec, group_inf, event_rate)])
        }
    }
    message(" - Final t = ", signif(epi_time, 5), ", values are:")
    print(table(pop$status))

    # tidy up pop
	pop[donor == 0 & status == "R", generation := 2]
    pop[, c("group_inf", "event_rate") := NULL]

    traits2 <- rbind(traits[sdp != "progeny"], pop, fill = TRUE)

    return(traits2)
}
