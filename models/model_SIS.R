#' Simulate an SIS model
#'
#' @param popn A population with trait values in data.table format
#' @param params A list of parameters
#'
#' @returns A new popn with epidemic event times, generation, and infectors
#' @export

# Main model ----
model_SIS <- function(popn, params) {
    message("Simulating an SIS epidemic ...")

    # copy necessary parameters
    r_beta <- params$r_beta
    t_end  <- params$t_end
    DEBUG  <- params$DEBUG

    # initialise population ----
    Xgroups <- init_popn(popn, params) |>
        split(by = "group")

    # function for appending to list unless only element is NA
    nax <- function(x, y) {
        if (length(x) == 1 && is.na(x[[1]])) y else c(x, y)
    }

    Y <- map(Xgroups, \(X) {
        if (FALSE) {
            X <- copy(Xgroups[[1]])
        }

        # This is unique to the SIS and SIRS models
        X[, `:=`(Tinf = as.list(Tinf),
                 Tdeath = as.list(Tdeath),
                 generation = as.list(generation),
                 infected_by = as.list(infected_by))]

        # This is a priority queue for the next event
        ni_events <- X[, .(.I, map_dbl(Tdeath, last))] |>
            melt(id.vars = "I", variable.name = "event", value.name = "time") |>
            setorder(time, na.last = TRUE)
        ni_events[, event := NULL]

        # Start epidemic simulation loop
        epi_time <- 0.0
        while (epi_time < t_end && X[, sum(status == "I") > 0L]) {
            if (DEBUG) message("time = ", signif(epi_time, 5))

            # Calculate infection rates in each group
            # this is the sum of the log infectivitities
            X[, group_inf := r_beta * GE * mean(inf * (status == "I"))]

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

            # check if next event is infection or non-infection
            if (epi_time + dt < t_next_event) {
                if (DEBUG) message("Next event is infection at t = ", signif(epi_time, 5))

                epi_time <- epi_time + dt

                # randomly select individual
                id_next_event <- sample(nrow(X), 1L, prob = X$inf_rate)

                xi <- X[, sample(.N, 1L, prob = inf * (status == "I"))]
                infd_by <- X$id[[xi]]
                next_gen <- X$generation[[xi]] + 1L

                sis_path <- generate_sis_path(epi_time, X, id_next_event, params)

                set(X, id_next_event,
                    c("status", "Tinf", "Tdeath", "generation", "infected_by"),
                    list("I",
                         nax(X$Tinf[[id_next_event]], epi_time),
                         nax(X$Tdeath[[id_next_event]], sis_path),
                         nax(X$generation[[id_next_event]], next_gen),
                         nax(X$infected_by[[id_next_event]], infd_by)))

                ni_events[I == id_next_event, time := sis_path]

                if (DEBUG) message("ID ", id_next_event, ": S -> I, infected by ID ", infd_by)
            } else {
                if (DEBUG) message("Next event is non-infection at t = ", signif(t_next_event, 5))

                epi_time <- t_next_event

                set(ni_events, 1L, "time", NA)

                status <- X[id_next_event, status]

                if (status == "I") {
                    if (DEBUG) message("ID ", id_next_event, ": I -> S")
                    set(X, id_next_event, "status", "S")
                } else {
                    print(X[, .(group, donor, status, Tinf, Tdeath, group_inf, inf_rate)])
                    print(X[id_next_event, .(group, donor, status, Tinf, Tdeath, group_inf, inf_rate)])
                    stop(str_c("selected ID ", id_next_event, "... unexpected event!"))
                    break
                }
            }

            # Every event we set the order of the ni_events. This should be fast
            # though as data.table uses a sensible sort method.
            setorder(ni_events, time, na.last = TRUE)

            if (DEBUG) {
                print(X[, .(group, donor, status, Tinf, Tdeath, group_inf, inf_rate)])
            }
        }
        X
    }) |> rbindlist()

    final_t <- Y$Tdeath |> list_c() |> max(na.rm = TRUE) |> signif(5)

    message(str_glue("- Final t = {final_t}, values are:"),
            str_flatten(capture.output(table(Y$status)), "\n"))

    # tidy up Y
    Y[, c("group_inf", "inf_rate") := NULL]
    Y[, parasites := !is.na(map_dbl(Tinf, last))]
    setorder(Y, id)

    popn2 <- rbind(popn[sdp != "progeny"], Y, fill = TRUE)

    return(popn2)
}


# A Susceptible individual's future disease trajectory is fixed at the point of
# exposure.
generate_sis_path <- function(epi_time, X, id, params) {
    with(params,
        rgamma(1L, RP_shape,
               scale = RP_scale * X$tol[[id]]) + epi_time)
}

