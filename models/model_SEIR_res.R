# Find the no. of individuals capable of infecting susceptibles at some point
# (so includes those currently exposed and not yet infectious).
get_seir_res_infectives <- function(X) {
    X[, sum(status %in% c("E", "I"))]
}


# A Susceptible individual's future disease trajectory is fixed at the point of
# exposure.
generate_seir_res_path <- function(epi_time, X, id, params) {
    with(params, {
        Tinf   <- epi_time
        Tsym   <- Tinf + rgamma(1L, LP_shape, scale = LP_scale * X$lat[[id]])
        Tdeath <- Tsym + rgamma(1L, DP_shape, scale = RP_scale * X$tol[[id]])

        # return
        list("E", Tinf, Tsym, Tdeath)
    })
}


# Find the time of the next non-infection event, and the id of the individual
next_seir_res_ni_event <- function(X, epi_time) {
    as.list(X[, .(.I,
                  Tsym2   = fifelse(Tsym > epi_time, Tsym, Inf),
                  Tdeath2 = fifelse(Tdeath > epi_time, Tdeath, Inf))]
            [, .(I, Tmin = pmin(Tsym2, Tdeath2))]
            [, list(t_next_event = min(Tmin, na.rm = TRUE), id_next_event = which.min(Tmin))])
}


# Main model ----
model_SEIR_res <- function(popn, params) {
    message("Simulating an SEIR_res epidemic ...")

    # copy necessary parameters
    r_beta   <- params$r_beta
    r_lambda <- 0.1 #params$r_lambda
    DEBUG    <- params$DEBUG

    # initialise populations ----
    X <- init_popn(popn, params)
    # group_inf is now the reservoir
    X[, group_inf := 0.0]

    # Start epidemic simulation loop
    dt <- 0.0
    epi_time <- 0.0
    while (get_seir_res_infectives(X) > 0L) {
        if (DEBUG) message("time = ", signif(epi_time, 5))

        # Calculate infection rates in each group
        # this is the sum of the log infectivitities
        X[, group_inf := group_inf * exp(- r_lambda * dt) +
                100 * GE * mean(fifelse(status == "I", inf, 0.0)) * r_beta * dt,
            by = group][]
        if (DEBUG) print(X[, mean(group_inf), by = group])

        # if S, infection at rate beta SI
        X[, event_rate := fifelse(status == "S", sus * group_inf, 0.0)]

        # id and time of next non-infection event
        ni_event      <- next_seir_res_ni_event(X, epi_time)
        t_next_event  <- ni_event$t_next_event
        id_next_event <- ni_event$id_next_event

        if (DEBUG) message("Next NI event id = ", id_next_event, " at t = ", t_next_event)

        # generate random timestep ----
        total_event_rate <- sum(X$event_rate)

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
            id_next_event <- sample(nrow(X),
                                    size = 1L,
                                    prob = X$event_rate)

            group_id <- X$group[id_next_event]
            infectives <- X[, .(.I, group, status, inf)][group == group_id & status == "I"]

            set(X, id_next_event, c("status", "Tinf", "Tsym", "Tdeath"),
                generate_seir_res_path(epi_time, X, id_next_event, params))

            if (DEBUG) message("ID ", id_next_event, ": S -> E")
        } else {
            if (DEBUG) message("next event is non-infection at t = ", signif(t_next_event, 5))

            dt <- t_next_event - epi_time
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
    # X[donor == 0 & status == "R", generation := 2]
    X[, c("group_inf", "event_rate") := NULL]
    Y[, parasites := !is.na(Tinf)]

    popn2 <- rbind(popn[sdp != "progeny"], X, fill = TRUE)

    popn2
}
