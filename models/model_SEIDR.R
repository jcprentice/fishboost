# Find the no. of individuals capable of infecting susceptibles at some point
# (so includes those currently exposed and not yet infectious).
get_seidr_infectives <- function(pop) {
    pop[, sum(status %in% c("E", "I", "D"))]
}


# A Susceptible individual's future disease trajectory is fixed at the point of
# exposure.
generate_seidr_path <- function(epi_time, pop, id, params) {
    with(params, {
        Tinf <- epi_time
        Tinc <- Tinf + rgamma(1L, r_eta_shape, r_eta_rate) * pop$latency[id]
        Tsym <- Tinc + rgamma(1L, r_rho_shape, r_rho_rate) * pop$detectability[id]
        Trec <- Tsym + rgamma(1L, r_gamma_shape, r_gamma_rate) * pop$recoverability[id]

        list("E", Tinf, Tinc, Tsym, Trec)
    })
}


# Find the time of the next non-infection event, and the id of the individual
next_seidr_ni_event <- function(pop, epi_time) {
    as.list(pop[, .(.I,
                    Tinc2 = fifelse(Tinc > epi_time, Tinc, Inf),
                    Tsym2 = fifelse(Tsym > epi_time, Tsym, Inf),
                    Trec2 = fifelse(Trec > epi_time, Trec, Inf))
    ][, .(I, Tmin = pmin(Tinc2, Tsym2, Trec2))
    ][, list(t_next_event = min(Tmin, na.rm = TRUE), id_next_event = which.min(Tmin))])
}


# Main model ----
model_SEIDR <- function(traits, params) {
    message("Simulating an SEIDR epidemic ...")

    # copy necessary parameters
    r_beta <- params$r_beta
    DEBUG  <- params$DEBUG

    # initialise populations ----
    pop <- init_pop(traits, params)

    # Start epidemic simulation loop ----
    epi_time <- 0.0
    while (get_seidr_infectives(pop) > 0L) {
        if (DEBUG) message("time = ", signif(epi_time, 5))

        # Calculate infection rates in each group ----
        # this is the sum of the log infectivitities
        # pop[, group_inf := mean(infectivity * r_beta * (status == "I")), by = group]
        pop[, group_inf := r_beta * GE * mean(fifelse(status %in% c("I", "D"), infectivity, 0.0)), by = group]

        # if S, infection at rate beta SI
        pop[, event_rate := fifelse(status == "S", susceptibility * group_inf, 0.0)]

        # id and time of next non-infection event
        ni_event      <- next_seidr_ni_event(pop, epi_time)
        t_next_event  <- ni_event$t_next_event
        id_next_event <- ni_event$id_next_event

        if (is.na(t_next_event)) t_next_event <- Inf

        if (DEBUG) message("next NI event id = ", id_next_event, " at t = ", signif(t_next_event, 5))

        # generate random timestep ----
        total_event_rate <- sum(pop$event_rate)

        # calculate dt if infections event rate > 0
        if (total_event_rate > 0.0) {
            dt <- rexp(1L, rate = total_event_rate)
        } else {
            dt <- Inf
        }

        if (DEBUG) message("Total infections event rate = ", signif(total_event_rate, 5))

        # check if next event is infection or non-infection ----
        if (epi_time + dt < t_next_event) {
            if (DEBUG) message("next event is infection at t = ", signif(epi_time, 5))

            epi_time <- epi_time + dt

            # randomly select individual
            id_next_event <- sample(nrow(pop),
                                    size = 1L,
                                    prob = pop$event_rate)

            group_id <- pop$group[id_next_event]
            infectives <- pop[, .(.I, group, status, infectivity)][group == group_id & status %in% c("I", "D")]
            infd_by <- safe_sample(x = infectives$I,
                                   size = 1L,
                                   prob = infectives$infectivity)
            next_gen <- pop$generation[infd_by] + 1L

            set(pop, id_next_event, c("status", "Tinf", "Tinc", "Tsym", "Trec"),
                generate_seidr_path(epi_time, pop, id_next_event, params))
            set(pop, id_next_event, c("generation", "infected_by"), list(next_gen, infd_by))

            if (DEBUG) message("ID ", id_next_event, ": S -> E, infected by ID ", infd_by)
        } else {
            if (DEBUG) message("next event is non-infection at t = ", signif(t_next_event, 5))

            epi_time <- t_next_event

            status <- pop$status[id_next_event]

            if (status == "E") {
                if (DEBUG) message("ID ", id_next_event, ": E -> I")
                set(pop, id_next_event, "status", "I")
            } else if (status == "I") {
                if (DEBUG) message("ID ", id_next_event, ": I -> D")
                set(pop, id_next_event, "status", "D")
            } else if (status == "D") {
                if (DEBUG) message("ID ", id_next_event, ": D -> R")
                set(pop, id_next_event, "status", "R")
            } else {
                message("status = ", status)
                print(pop[, .(group, donor, status, Tinf, Tinc, Tsym, Trec, group_inf, event_rate)])
                print(pop[id_next_event, .(group, donor, status, Tinf, Tinc, Tsym, Trec, group_inf, event_rate)])
                stop("selected ID ", id_next_event, "... unexpected event!")
                break
            }
        }

        if (DEBUG) {
            print(pop[, .(group, donor, status, Tinf, Tinc, Tsym, Trec, group_inf, event_rate)])
        }
    }
    message(" - Final t = ", signif(epi_time, 5), ", values are:")
    print(table(pop$status))

    # tidy up pop
    pop[, c("group_inf", "event_rate") := NULL]

    # fix generation
    pop[status == "R" & donor == 0L, generation := 2L]

    traits2 <- rbind(traits[sdp != "progeny"], pop, fill = TRUE)

    return(traits2)
}
