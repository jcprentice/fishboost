# Find the no. of individuals capable of infecting susceptibles at some point
# (so includes those currently exposed and not yet infectious).
get_seidr_infectives <- function(X) {
    X[, sum(status %in% c("E", "I", "D"))]
}


# A Susceptible individual's future disease trajectory is fixed at the point of
# exposure.
generate_seidr_path <- function(epi_time, X, id, params) {
    with(params, {
        Tinf   <- epi_time
        Tinc   <- Tinf + rgamma(1L, LP_shape, scale = LP_scale * X$lat[[id]])
        Tsym   <- Tinc + rgamma(1L, DP_shape, scale = DP_scale * X$det[[id]])
        Tdeath <- Tsym + rgamma(1L, RP_shape, scale = RP_scale * X$tol[[id]])

        list("E", Tinf, Tinc, Tsym, Tdeath)
    })
}


# Find the time of the next non-infection event, and the id of the individual
next_seidr_ni_event <- function(X, epi_time) {
    as.list(X[, .(.I,
                  Tinc2 = fifelse(Tinc > epi_time, Tinc, Inf),
                  Tsym2 = fifelse(Tsym > epi_time, Tsym, Inf),
                  Tdeath2 = fifelse(Tdeath > epi_time, Tdeath, Inf))]
            [, .(I, Tmin = pmin(Tinc2, Tsym2, Tdeath2))]
            [, list(t_next_event = min(Tmin, na.rm = TRUE), id_next_event = which.min(Tmin))])
}


# Main model ----
model_SEIDR <- function(popn, params) {
    message("Simulating an SEIDR epidemic ...")

    # copy necessary parameters
    r_beta <- params$r_beta
    DEBUG  <- params$DEBUG

    # initialise populations ----
    X <- init_popn(popn, params)

    # Start epidemic simulation loop ----
    epi_time <- 0.0
    while (get_seidr_infectives(X) > 0L) {
        if (DEBUG) message("time = ", signif(epi_time, 5))

        # Calculate infection rates in each group ----
        # this is the sum of the log infectivitities
        # X[, group_inf := mean(inf * r_beta * (status == "I")), by = group]
        X[, group_inf := r_beta * GE * mean(fifelse(status %in% c("I", "D"), inf, 0.0)), by = group]
        
        if (DEBUG) {
            message("Group infections are:\n",
                    str_flatten(capture.output(X[, .(inf = signif(group_inf[1], 3)), group]), "\n"))
        }

        # if S, infection at rate beta SI
        X[, event_rate := fifelse(status == "S", sus * group_inf, 0.0)]

        # id and time of next non-infection event
        ni_event      <- next_seidr_ni_event(X, epi_time)
        t_next_event  <- ni_event$t_next_event
        id_next_event <- ni_event$id_next_event

        if (is.na(t_next_event)) t_next_event <- Inf

        if (DEBUG) message("Next NI event id = ", id_next_event, " at t = ", signif(t_next_event, 5))

        # generate random timestep ----
        total_event_rate <- sum(X$event_rate)

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
            id_next_event <- sample(nrow(X),
                                    size = 1L,
                                    prob = X$event_rate)

            group_id <- X$group[id_next_event]
            infectives <- X[, .(.I, group, status, inf)][group == group_id & status %in% c("I", "D")]
            infd_by <- safe_sample(x = infectives$I,
                                   size = 1L,
                                   prob = infectives$inf)
            next_gen <- X$generation[infd_by] + 1L
            
            if (DEBUG) {
                message("group infectivities = \n",
                        str_flatten(capture.output(X[, .(inf = group_inf[1]), group]), "\n"))
            }

            set(X, id_next_event, c("status", "Tinf", "Tinc", "Tsym", "Tdeath"),
                generate_seidr_path(epi_time, X, id_next_event, params))
            set(X, id_next_event, c("generation", "infected_by"), list(next_gen, infd_by))

            if (DEBUG) message("ID ", id_next_event, ": S -> E, infected by ID ", infd_by)
        } else {
            if (DEBUG) message("next event is non-infection at t = ", signif(t_next_event, 5))

            epi_time <- t_next_event

            status <- X$status[id_next_event]

            if (status == "E") {
                if (DEBUG) message("ID ", id_next_event, ": E -> I")
                set(X, id_next_event, "status", "I")
            } else if (status == "I") {
                if (DEBUG) message("ID ", id_next_event, ": I -> D")
                set(X, id_next_event, "status", "D")
            } else if (status == "D") {
                if (DEBUG) message("ID ", id_next_event, ": D -> R")
                set(X, id_next_event, "status", "R")
            } else {
                message("status = ", status)
                print(X[, .(group, donor, status, Tinf, Tinc, Tsym, Tdeath, group_inf, event_rate)])
                print(X[id_next_event, .(group, donor, status, Tinf, Tinc, Tsym, Tdeath, group_inf, event_rate)])
                stop("selected ID ", id_next_event, "... unexpected event!")
                break
            }
        }

        if (DEBUG == 2) {
            print(X[, .(group, donor, status, Tinf, Tinc, Tsym, Tdeath, group_inf, event_rate)])
        }
    }
    message(" - Final t = ", signif(epi_time, 5), ", values are:",
            str_flatten(capture.output(table(X$status)), "\n"))

    # tidy up X
    X[, c("group_inf", "event_rate") := NULL]

    # fix generation
    # X[status == "R" & donor == 0L, generation := 2L]

    popn2 <- rbind(popn[sdp != "progeny"], X, fill = TRUE)

    popn2
}
