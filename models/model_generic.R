# Find the no. of individuals capable of infecting susceptibles at some point
# (so includes those currently exposed and not yet infectious).
get_infectives <- function(X, params) {
    X[, sum(status %in% params$infectious_subset)]
}


# A Susceptible individual's future disease trajectory is fixed at the point of
# exposure.
generate_path <- function(epi_time, X, id, params) {
    with(params, {
        Tinf <- epi_time

        if (model_type %in% c("SEIDR", "SEIDR_res")) {
            Tinc   <- Tinf + rgamma(1L, LP_shape, scale = LP_scale * X$lat[[id]])
            Tsym   <- Tinc + rgamma(1L, DP_shape, scale = DP_scale * X$det[[id]])
            Tdeath <- Tsym + rgamma(1L, RP_shape, scale = RP_scale * X$tol[[id]])
            list("E", Tinf, Tinc, Tsym, Tdeath)
        } else if (model_type %in% c("SEIR", "SEIR_res")) {
            Tsym   <- Tinf + rgamma(1L, LP_shape, scale = LP_scale * X$lat[[id]])
            Tdeath <- Tsym + rgamma(1L, RP_shape, scale = RP_scale * X$tol[[id]])
            list("E", Tinf, Tsym, Tdeath)
        } else if (model_type %in% c("SIDR", "SIDR_res")) {
            Tsym   <- Tinf + rgamma(1L, DP_shape, scale = DP_scale * X$det[[id]])
            Tdeath <- Tsym + rgamma(1L, RP_shape, scale = RP_scale * X$tol[[id]])
            list("I", Tinf, Tsym, Tdeath)
        } else if (model_type %in% c("SIS", "SIR", "SIS_res", "SIR_res")) {
            Tdeath <- Tinf + rgamma(1L, RP_shape, scale = RP_scale * X$tol[[id]])
            list("I", Tinf, Tdeath)
        }
    })
}


# Find the time of the next non-infection event, and the id of the individual
next_ni_event <- function(X, epi_time) {
    as.list(X[, .(.I,
                  Tsym2   = fifelse(Tsym   > epi_time, Tsym,   Inf),
                  Tdeath2 = fifelse(Tdeath > epi_time, Tdeath, Inf))]
            [, .(I, Tmin = pmin(Tsym2, Tdeath2))]
            [, list(t_next_event = min(Tmin, na.rm = TRUE), id_next_event = which.min(Tmin))])
}


get_infection_rates <- function(X, dt, params) {
    r_beta <- params$r_beta

    X[, group_inf := r_beta * GE * mean(fifelse(status == "I", inf, 0.0)), by = group]

    if (params$reservoir) {
        X[, group_res := group_res * exp(- params$r_lambda * dt) +
                100 * GE * mean(fifelse(status == "I", inf, 0.0)) * r_beta * dt,
            by = group]
    }
}


# Main model ----
model_generic <- function(popn, params) {

    # copy necessary parameters
    model_type <- params$model_type
    r_beta     <- params$r_beta
    DEBUG      <- params$DEBUG


    message("Simulating an ", model_type, " epidemic ...")

    # initialise populations ----
    X <- init_popn(popn, params)

    # Start epidemic simulation loop ----
    epi_time <- dt <- 0.0
    while (get_infectives(X) > 0L) {
        if (DEBUG) message("time = ", signif(epi_time, 5L))

        # Calculate infection rates in each group ----
        # this is the sum of the log infectivitities
        get_infection_rates(X, params)

        if (DEBUG) {
            message("Group infections are:\n",
                    str_flatten(capture.output(X[, signif(group_inf[1], 3), group]), "\n"))
        }

        # if S, infection at rate beta SI
        X[, event_rate := fifelse(status == "S", sus * (group_inf + group_res), 0.0)]

        # id and time of next non-infection event
        ni_event      <- next_ni_event(X, epi_time)
        t_next_event  <- ni_event$t_next_event
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
            infectives <- X[, .(.I, group, status, inf)][group == group_id & status == "I"]
            infd_by <- safe_sample(x = infectives$I,
                                   size = 1L,
                                   prob = infectives$inf)
            next_gen <- X[infd_by, generation + 1L]

            set(X, id_next_event, c("status", "Tinf", "Tsym", "Tdeath"),
                generate_path(epi_time, X, id_next_event, params))
            set(X, id_next_event, c("generation", "infected_by"), list(next_gen, infd_by))

            if (DEBUG) message("ID ", id_next_event, ": S -> E, infected by ID ", infd_by)
        } else {
            if (DEBUG) message("next event is non-infection at t = ", signif(t_next_event, 5))

            epi_time <- t_next_event

            status <- X[id_next_event, status]

            if (status == "E") {
                if (DEBUG) message("ID ", id_next_event, ": E -> I")
                set(X, id_next_event, "status", "I")
            } else if (status == "I") {
                if (DEBUG) message("ID ", id_next_event, ": I -> R")
                set(X, id_next_event, "status", "R")
            } else {
                print(X[, .(group, donor, status, Tinf, Tsym, Tdeath, group_inf, event_rate)])
                print(X[id_next_event, .(group, donor, status, Tinf, Tsym, Tdeath, group_inf, event_rate)])
                stop(str_c("selected ID ", id_next_event, "... unexpected event!"))
                break
            }
        }

        if (DEBUG) {
            print(X[, .(group, donor, status, Tinf, Tsym, Tdeath, group_inf, event_rate)])
        }
    }

    final_t <- epi_time |> signif(5)

    message(str_glue("- Final t = {final_t}, values are:"),
            str_flatten(capture.output(table(X$status)), "\n"))

    # tidy up X
    X[, c("group_inf", "group_res", "event_rate") := NULL]
    X[, parasites := !is.na(Tinf)]

    # fix generation
    # X[status == "R" & donor == 0, generation := 2]

    popn2 <- rbind(popn[sdp != "progeny"], X, fill = TRUE)

    return(popn2)
}
