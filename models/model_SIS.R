# Find the no. of individuals capable of infecting susceptibles at some point
# (so includes those currently exposed and not yet infectious).
get_sis_infectives <- function(X) {
    X[, sum(status == "I")]
}


# A Susceptible individual's future disease trajectory is fixed at the point of
# exposure.
generate_sis_path <- function(epi_time, X, id, params) {
    with(params, {
        Tinf   <- epi_time
        Tdeath <- Tinf + rgamma(1L, RP_shape, scale = RP_scale * X$tol[[id]])

        list("I", Tinf, Tdeath)
    })
}


# Find the time of the next non-infection event, and the id of the individual
next_sis_ni_event <- function(X, epi_time) {
    as.list(X[, .(.I,
                  Tmin = fifelse(Tdeath > epi_time, Tdeath, Inf))]
            [, list(t_next_event = min(Tmin, na.rm = TRUE), id_next_event = which.min(Tmin))])
}


# Main model ----
model_SIS <- function(popn, params) {
    message("Simulating an SIS epidemic ...")

    r_beta <- params$r_beta
    tmax   <- params$tmax
    DEBUG  <- params$DEBUG

    X <- copy(popn)

    # initialise population ----
    X <- init_popn(popn, params)

    # Start epidemic simulation loop
    epi_time <- 0.0
    while (get_sis_infectives(X) > 0L && epi_time < tmax) {
        if (DEBUG) message("time = ", signif(epi_time, 5))

        # Calculate infection rates in each group
        # this is the sum of the log infectivitities
        X[, group_inf := r_beta * GE * mean(inf * (status == "I")), by = group]

        # if S, infection at rate beta SI
        X[, event_rate := fifelse(status == "S", sus * group_inf, 0.0)]

        # id and time of next non-infection event
        ni_event <- next_sis_ni_event(X, epi_time)
        t_next_event <- ni_event$t_next_event
        id_next_event <- ni_event$id_next_event

        if (DEBUG) message("Next NI event id = ", id_next_event, " at t = ", t_next_event)


        # generate random timestep ----
        total_event_rate <- sum(X$event_rate)

        # calculate dt if infections event rate > 0
        if (total_event_rate > 0.0) {
            dt <- rexp(1L, rate = total_event_rate)
        } else {
            dt <- Inf
        }
        if (DEBUG) message("Total Infections Event Rate = ", signif(total_event_rate, 5))

        # check if next event is infection or non-infection
        if (epi_time + dt < t_next_event) {
            if (DEBUG) message("next event is infection at t = ", signif(epi_time, 5))

            epi_time <- epi_time + dt

            # randomly select individual
            id_next_event <- sample(nrow(X),
                                    size = 1L,
                                    prob = X$event_rate)

            group_id <- X$group[id_next_event]
            infectives <- X[, .(.I, group, status, inf)][group == group_id & status == "I"]
            infd_by <- safe_sample(x = infectives$I,
                                   size = 1L,
                                   prob = infectives$inf)
            next_gen <- X[infd_by, generation + 1L]

            set(X, id_next_event, c("status", "Tinf", "Tdeath"),
                generate_sir_path(epi_time, X, id_next_event, params))
            set(X, id_next_event, c("generation", "infected_by"), list(next_gen, infd_by))

            if (DEBUG) message("ID ", id_next_event, ": S -> I, infected by ID ", infd_by)
        } else {
            if (DEBUG) message("next event is non-infection at t = ", signif(t_next_event, 5))

            epi_time <- t_next_event

            status <- X[id_next_event, status]

            if (status == "I") {
                if (DEBUG) message("ID ", id_next_event, ": I -> S")
                set(X, id_next_event, "status", "S")
            } else {
                print(X[, .(group, donor, status, Tinf, Tdeath, group_inf, event_rate)])
                print(X[id_next_event, .(group, donor, status, Tinf, Tdeath, group_inf, event_rate)])
                stop(str_c("selected ID ", id_next_event, "... unexpected event!"))
                break
            }
        }

        if (DEBUG) {
            print(X[, .(group, donor, status, Tinf, Tdeath, group_inf, event_rate)])
        }
    }
    message(" - Final t = ", signif(epi_time, 5), ", values are:",
            str_flatten(capture.output(table(X$status)), "\n"))

    # tidy up X
    X[, c("group_inf", "event_rate") := NULL]
    X[, parasites := !is.na(Tinf)]
    setorder(X, id)

    popn2 <- rbind(popn[sdp != "progeny"], X, fill = TRUE)

    return(popn2)
}
