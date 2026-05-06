#' Simulate an SIDR model
#'
#' @param popn A population with trait values in data.table format
#' @param params A list of parameters
#'
#' @returns A new popn with epidemic event times, generation, and infectors
#' @export

# Main model ----
model_SIDR <- function(popn, params) {
    message("Simulating an SIDR groups epidemic ...")

    # copy necessary parameters
    r_beta <- params$r_beta
    DEBUG  <- params$DEBUG

    # initialise populations ----
    Xgroups <- init_popn(popn, params) |>
        split(by = "group")

    # Each group is independent, so it seems to be about ~33% faster to split
    # them by group, run them separately, then recombine them.

    Y <- map(Xgroups, \(X) {
        # This is a priority queue for the next event
        ni_events <- X[, .(.I, Tsym, Tdeath)] |>
            melt(id.vars = "I", variable.name = "event", value.name = "time") |>
            setorder(time, na.last = TRUE)
        ni_events[, event := NULL]

        # Start epidemic simulation loop ----
        epi_time <- 0.0
        while (X[, sum(status %in% c("I", "D")) > 0]) {
            if (DEBUG) message("time = ", signif(epi_time, 5))

            # Calculate infection rates in each group ----
            X[, group_inf := r_beta * GE * mean(fifelse(status %in% c("I", "D"), inf, 0.0))]

            # if S, infection at rate beta SI
            X[, inf_rate := fifelse(status == "S", sus * group_inf, 0.0)]

            # id and time of next non-infection event
            t_next_event  <- ni_events$time[[1]]
            id_next_event <- ni_events$I[[1]]

            if (is.na(t_next_event)) t_next_event <- Inf

            if (DEBUG) message("Next NI event id = ", id_next_event, " at t = ", signif(t_next_event, 5))

            # generate random timestep ----
            total_inf_rate <- sum(X$inf_rate)

            # calculate dt if infections event rate > 0
            dt <- if (total_inf_rate > 0.0) {
                rexp(1L, rate = total_inf_rate)
            } else {
                Inf
            }

            if (DEBUG) message("Total infections event rate = ", signif(total_inf_rate, 5))

            # check if next event is infection or non-infection ----
            if (epi_time + dt < t_next_event) {
                if (DEBUG) message("Next event is infection at t = ", signif(epi_time, 5))

                epi_time <- epi_time + dt

                # randomly select individual
                id_next_event <- samp1(X$inf_rate)

                infectives <- X[status %in% c("I", "D"), .(.I, status, inf)]
                infd_by <- safe_sample(x = infectives$I,
                                       size = 1L,
                                       prob = infectives$inf)
                next_gen <- X$generation[[infd_by]] + 1L

                sidr_path <- generate_sidr_path(epi_time, X, id_next_event, params)
                set(X, id_next_event, c("status", "Tinf", "Tsym", "Tdeath"), sidr_path)

                ni_events[I == id_next_event, time := as.numeric(sidr_path[-(1:2)])]

                set(X, id_next_event, c("generation", "infected_by"), list(next_gen, infd_by))

                if (DEBUG) message("ID ", id_next_event, ": S -> I, infected by ID ", infd_by)
            } else {
                if (DEBUG) message("Next event is non-infection at t = ", signif(t_next_event, 5))

                epi_time <- t_next_event

                set(ni_events, 1L, "time", NA)

                status <- X$status[[id_next_event]]

                if (status == "I") {
                    if (DEBUG) message("ID ", id_next_event, ": I -> D")
                    set(X, id_next_event, "status", "D")
                } else if (status == "D") {
                    if (DEBUG) message("ID ", id_next_event, ": D -> R")
                    set(X, id_next_event, "status", "R")
                } else {
                    message("status = ", status)
                    print(X[, .(group, donor, status, Tinf, Tsym, Tdeath, group_inf, inf_rate)])
                    stop("selected ID ", id_next_event, "... unexpected event!")
                    break
                }
            }

            # Every event we set the order of the ni_events. This should be fast
            # though as data.table uses a sensible sort method.
            setorder(ni_events, time, na.last = TRUE)

            if (DEBUG == 2) {
                print(X[, .(group, donor, status, Tinf, Tsym, Tdeath, group_inf, inf_rate)])
            }
        }
        X
    }) |> rbindlist()

    final_t <- max(Y$Tdeath, na.rm = TRUE) |> signif(5)

    message(str_glue("- Final t = {final_t}, values are:"),
            str_flatten(capture.output(table(Y$status)), "\n"))

    # tidy up X
    Y[, c("group_inf", "inf_rate") := NULL]
    Y[, parasites := !is.na(Tinf)]

    # fix generation
    # X[status == "R" & donor == 0L, generation := 2L]

    popn2 <- rbind(popn[sdp != "progeny"], Y, fill = TRUE)

    popn2
}

# A Susceptible individual's future disease trajectory is fixed at the point of
# exposure.
generate_sidr_path <- function(epi_time, X, id, params) {
    with(params, {
        Tinf   <- epi_time
        Tsym   <- Tinf + rgamma(1L, DP_shape, scale = DP_scale * X$det[[id]])
        Tdeath <- Tsym + rgamma(1L, RP_shape, scale = RP_scale * X$tol[[id]])

        list("I", Tinf, Tsym, Tdeath)
    })
}
