# Find the no. of individuals capable of infecting susceptibles at some point
# (so includes those currently exposed and not yet infectious).
get_sir_infectives <- function(X) {
    X[, sum(status == "I")]
}


# A Susceptible individual's future disease trajectory is fixed at the point of
# exposure.
generate_sir_path <- function(epi_time, X, id, params) {
    with(params, {
        Tinf   <- epi_time
        Tdeath <- Tinf + rgamma(1L, RP_shape, scale = RP_scale * X$tol[[id]])

        list("I", Tinf, Tdeath)
    })
}


# Find the time of the next non-infection event, and the id of the individual
next_sir_ni_event <- function(X, epi_time) {
    as.list(X[, .(.I,
                  Tmin = fifelse(Tdeath > epi_time, Tdeath, Inf))]
            [, list(t_next_event = min(Tmin, na.rm = TRUE), id_next_event = which.min(Tmin))])
}


# Main model ----
model_SIR <- function(popn, params) {
    message("Simulating an SIR epidemic ...")

    # copy necessary parameters
    r_beta <- params$r_beta
    DEBUG  <- params$DEBUG

    # initialise population ----
    X <- init_popn(popn, params)

    # Start epidemic simulation loop ----
    epi_time <- 0.0
    while (get_sir_infectives(popn) > 0L) {
        if (DEBUG) message("time = ", signif(epi_time, 5))

        # Calculate infection rates in each group ----
        # this is the sum of the log infectivitities
        popn[, group_inf := r_beta * GE * mean(inf * (status == "I")), by = group]

        # if S, infection at rate beta SI
        popn[, event_rate := fifelse(status == "S", sus * group_inf, 0.0)]

        # id and time of next non-infection event
        ni_event <- next_sir_ni_event(popn, epi_time)
        t_next_event <- ni_event$t_next_event
        id_next_event <- ni_event$id_next_event

        if (DEBUG) message("Next NI event id = ", id_next_event, " at t = ", t_next_event)


        # generate random timestep ----
        total_event_rate <- sum(popn$event_rate)

        # calculate dt if infections event rate > 0
        if (total_event_rate > 0.0) {
            dt <- rexp(1L, rate = total_event_rate)
        } else {
            dt <- Inf
        }
        if (DEBUG) message("Total Infections Event Rate = ", signif(total_event_rate, 5))

        # check if next event is infection or non-infection ----
        if (epi_time + dt < t_next_event) {
            if (DEBUG) message("next event is infection at t = ", signif(epi_time, 5))

            epi_time <- epi_time + dt

            # randomly select individual
            id_next_event <- sample(nrow(popn),
                                    size = 1L,
                                    prob = popn$event_rate)

            group_id <- popn$group[id_next_event]
            # infectives <- popn[, .(.I, group, status, inf)][group == group_id & status == "I"]
            infectives <- popn[group == group_id & status == "I", .(.I, group, status, inf)]
            infd_by <- safe_sample(x = infectives$I,
                                   size = 1L,
                                   prob = infectives$inf)
            next_gen <- popn[infd_by, generation + 1L]

            set(popn, id_next_event, c("status", "Tinf", "Tdeath"),
                generate_sir_path(epi_time, popn, id_next_event, params))
            set(popn, id_next_event, c("generation", "infected_by"), list(next_gen, infd_by))

            if (DEBUG) message("ID ", id_next_event, ": S -> I, infected by ID ", infd_by)
        } else {
            if (DEBUG) message("next event is non-infection at t = ", signif(t_next_event, 5))

            epi_time <- t_next_event

            status <- popn[id_next_event, status]

            if (status == "I") {
                if (DEBUG) message("ID ", id_next_event, ": I -> R")
                set(popn, id_next_event, "status", "R")
            } else {
                print(popn[, .(group, donor, status, Tinf, Tdeath, group_inf, event_rate)])
                print(popn[id_next_event, .(group, donor, status, Tinf, Tdeath, group_inf, event_rate)])
                stop(str_c("selected ID ", id_next_event, "... unexpected event!"))
                break
            }
        }

        if (DEBUG) {
            print(popn[, .(group, donor, status, Tinf, Tdeath, group_inf, event_rate)])
        }
    }
    message(" - Final t = ", signif(epi_time, 5), ", values are:",
            str_flatten(capture.output(table(popn$status)), "\n"))

    # tidy up popn
    popn[, c("group_inf", "event_rate") := NULL]

    popn2 <- rbind(popn[sdp != "progeny"], X, fill = TRUE)

    popn2
}
