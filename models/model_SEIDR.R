#' Simulate an SEIDR model
#'
#' @param popn A population with trait values in data.table format
#' @param params A list of parameters
#'
#' @returns A new popn with epidemic event times, generation, and infectors
#' @export

# Main model ----
model_SEIDR <- function(popn, params) {
    message("Simulating an SEIDR groups epidemic ...")

    # copy necessary parameters
    r_beta <- params$r_beta
    DEBUG  <- params$DEBUG

    # initialise populations ----
    Xgroups <- init_popn(popn, params) |>
        split(by = "group")

    # Each group is independent, so it seems to be about ~33% faster to split
    # them by group, run them separately, then recombine them.

    Y <- map(Xgroups, \(X) {
        if (FALSE) {
            X <- copy(Xgroups[[1]])
        }
        # This is a priority queue for the next event
        ni_events <- X[, .(.I, Tinc, Tsign, Tdeath)] |>
            melt(id.vars = "I", variable.name = "event", value.name = "time") |>
            setorder(time, na.last = TRUE)
        ni_events[, event := NULL]

        # Start epidemic simulation loop ----
        epi_time <- 0.0
        while (X[, sum(status %in% c("E", "I", "D")) > 0]) {
            if (DEBUG) message("time = ", signif(epi_time, 5))

            # Calculate infection rates in each group ----
            X[, group_inf := r_beta * GE * mean(inf * (status %in% c("I", "D")))]

            # if S, infection at rate beta SI
            X[, inf_rate := sus * group_inf * (status == "S")]

            # id and time of next non-infection event
            t_next_event  <- ni_events$time[[1]]
            id_next_event <- ni_events$I[[1]]

            if (is.na(t_next_event)) t_next_event <- Inf

            if (DEBUG) message("Next NI event id = ", id_next_event,
                               " at t = ", signif(t_next_event, 5))

            # generate random timestep ----
            total_inf_rate <- sum(X$inf_rate)

            # calculate dt if infections event rate > 0
            dt <- rexp(1L) / total_inf_rate

            if (DEBUG) message("Total infections event rate = ", signif(total_inf_rate, 5))

            # check if next event is infection or non-infection ----
            if (epi_time + dt < t_next_event) {
                if (DEBUG) message("Next event is infection at t = ", signif(epi_time, 5))

                epi_time <- epi_time + dt

                # randomly select individual
                id_next_event <- sample(nrow(X), 1L, prob = X$inf_rate)

                xi <- X[, sample(.N, 1L, prob = inf * (status %in% c("I", "D")))]
                infd_by <- X$id[[xi]]
                next_gen <- X$generation[[xi]] + 1L

                seidr_path <- generate_seidr_path(epi_time, X, id_next_event, params)

                set(X, id_next_event,
                    c("status", "Tinf", "Tinc", "Tsign", "Tdeath", "generation", "infected_by"),
                    c("E", epi_time, as.list(seidr_path), next_gen, infd_by))

                ni_events[I == id_next_event, time := seidr_path]

                if (DEBUG) message("ID ", id_next_event, ": S -> E, infected by ID ", infd_by)
            } else {
                if (DEBUG) message("Next event is non-infection at t = ", signif(t_next_event, 5))

                epi_time <- t_next_event

                set(ni_events, 1L, "time", NA)

                status <- X[id_next_event, as.character(status)]

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
                    print(X[, .(group, donor, status, Tinf, Tinc, Tsign, Tdeath, group_inf, inf_rate)])
                    print(X[id_next_event, .(group, donor, status, Tinf, Tinc, Tsign, Tdeath, group_inf, inf_rate)])
                    stop("selected ID ", id_next_event, "... unexpected event!")
                    break
                }
            }

            # Every event we set the order of the ni_events. This should be fast
            # though as data.table uses a sensible sort method.
            setorder(ni_events, time, na.last = TRUE)

            if (DEBUG == 2) {
                print(X[, .(donor, status, Tinf, Tinc, Tsign, Tdeath,
                            generation, infected_by)])
                # print(X[, .(donor, status, generation, infected_by)])
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
generate_seidr_path <- function(epi_time, X, id, params) {
    with(params,
         rgamma(3L, c(LP_shape, DP_shape, RP_shape),
                scale = c(LP_scale * X$lat[[id]],
                          DP_scale * X$det[[id]],
                          RP_scale * X$tol[[id]])) |>
             cumsum() + epi_time)
}

