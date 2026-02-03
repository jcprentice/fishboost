make_time_series <- function(popn, params) {
    switch(params$model_type,
           "SIR" = make_time_series_sir(popn, params),
           "SEIR" = make_time_series_seir(popn, params),
           "SIDR" = make_time_series_sidr(popn, params),
           "SEIDR" = make_time_series_seidr(popn, params)
    )
}


make_time_series_seidr <- function(popn, params) {
    compartments <- c("S", "E", "I", "D", "R")
    # tmax <- params$tmax
    tmax <- max(popn$Tdeath, na.rm = TRUE)

    popn2 <- popn[sdp == "progeny", .(Tinf, Tinc, Tsym, Tdeath)]
    popn2[is.na(Tsym), Tsym := Tdeath]
    popn2[is.na(Tinc), Tinc := Tsym]
    popn2[is.na(Tinf), Tinf := Tinc]

    N <- popn2[, .N]

    X <- rbind(
        data.table(event = factor("start", c("start", "infection", "incubation", "detection", "removal", "end")),
                   time = 0.0),
        data.table(event = "infection",  time = popn2$Tinf[popn2$Tinf > 0]),
        data.table(event = "incubation", time = popn2$Tinc),
        data.table(event = "detection",  time = popn2$Tsym),
        data.table(event = "removal",    time = popn2$Tdeath),
        data.table(event = "end", time = tmax))

    X <- X[is.finite(time)]

    setorder(X, time)
    #                                                S    E    I    D    R
    X[event == "start",      (compartments) := list(+N,  +0L, +0L, +0L, +0L)]
    X[event == "infection",  (compartments) := list(-1L, +1L, +0L, +0L, +0L)]
    X[event == "incubation", (compartments) := list(+0L, -1L, +1L, +0L, +0L)]
    X[event == "detection",  (compartments) := list(+0L, +0L, -1L, +1L, +0L)]
    X[event == "removal",    (compartments) := list(+0L, +0L, +0L, -1L, +1L)]
    X[event == "end",        (compartments) := list(+0L, +0L, +0L, +0L, +0L)]

    X[, event := NULL]
    X[, S := cumsum(S)]
    X[, E := cumsum(E)]
    X[, I := cumsum(I)]
    X[, D := cumsum(D)]
    X[, R := cumsum(R)]
    X[, ID := I + D]

    X <- X[time <= tmax]

    X
}


make_time_series_sidr <- function(popn, params) {
    compartments <- c("S", "I", "D", "R")
    # tmax <- params$tmax
    tmax <- max(popn$Tdeath, na.rm = TRUE)

    popn2 <- popn[sdp == "progeny", .(Tinf, Tsym, Tdeath)]
    popn2[is.na(Tsym), Tsym := Tdeath]
    popn2[is.na(Tinf), Tinf := Tsym]

    N <- popn2[, .N]

    X <- rbind(
        data.table(event = factor("start", c("start", "infection", "detection", "removal", "end")),
                   time = 0.0),
        data.table(event = "infection", time = popn2$Tinf[popn2$Tinf > 0]),
        data.table(event = "detection", time = popn2$Tsym),
        data.table(event = "removal",   time = popn2$Tdeath),
        data.table(event = "end",       time = tmax))

    X <- X[is.finite(time)]

    setorder(X, time)
    #                                                S    I    D    R
    X[event == "start",     (compartments) := list(+N,  +0L, +0L, +0L)]
    X[event == "infection", (compartments) := list(-1L, +1L, +0L, +0L)]
    X[event == "detection", (compartments) := list(+0L, -1L, +1L, +0L)]
    X[event == "removal",   (compartments) := list(+0L, +0L, -1L, +1L)]
    X[event == "end",       (compartments) := list(+0L, +0L, +0L, +0L)]

    X[, event := NULL]
    X[, S := cumsum(S)]
    X[, I := cumsum(I)]
    X[, D := cumsum(D)]
    X[, R := cumsum(R)]
    X[, ID := I + D]

    X <- X[time <= tmax]

    X
}


make_time_series_seir <- function(popn, params) {
    compartments <- c("S", "E", "I", "R")
    # tmax <- params$tmax
    tmax <- max(popn$Tdeath, na.rm = TRUE)

    popn2 <- popn[sdp == "progeny", .(Tinf, Tsym, Tdeath)]
    popn2[is.na(Tsym), Tsym := Tdeath]
    popn2[is.na(Tinf), Tinf := Tsym]

    N <- popn2[, .N]

    X <- rbind(
        data.table(event = factor("start", c("start", "infection", "incubation", "removal", "end")),
                   time = 0.0),
        data.table(event = "infection",  time = popn2$Tinf[popn2$Tinf > 0]),
        data.table(event = "incubation", time = popn2$Tsym),
        data.table(event = "removal",    time = popn2$Tdeath),
        data.table(event = "end",        time = tmax))

    X <- X[is.finite(time)]

    setorder(X, time)

    X[event == "start",      (compartments) := list(+N,  +0L, +0L, +0L)]
    X[event == "infection",  (compartments) := list(-1L, +1L, +0L, +0L)]
    X[event == "incubation", (compartments) := list(+0L, -1L, +1L, +0L)]
    X[event == "removal",    (compartments) := list(+0L, +0L, -1L, +1L)]
    X[event == "end",        (compartments) := list(+0L, +0L, +0L, +0L)]

    X[, event := NULL]
    X[, S := cumsum(S)]
    X[, E := cumsum(E)]
    X[, I := cumsum(I)]
    X[, R := cumsum(R)]

    X <- X[time <= tmax]

    X
}


make_time_series_sir <- function(popn, params) {
    compartments <- c("S", "I", "R")
    # tmax <- params$tmax
    tmax <- max(popn$Tdeath, na.rm = TRUE)

    popn2 <- popn[sdp == "progeny", .(Tinf, Tdeath)]
    popn2[is.na(Tinf), Tinf := Tdeath]

    N <- popn2[, .N]

    X <- rbind(
        data.table(event = factor("start", c("start", "infection", "removal", "end")),
                   time = 0.0),
        data.table(event = "infection", time = popn2$Tinf[popn2$Tinf > 0]),
        data.table(event = "removal",   time = popn2$Tdeath),
        data.table(event = "end",       time = tmax))

    X <- X[is.finite(time)]

    setorder(X, time)

    X[event == "start",     (compartments) := list(+N,  +0L, +0L)]
    X[event == "infection", (compartments) := list(-1L, +1L, +0L)]
    X[event == "removal",   (compartments) := list(+0L, -1L, +1L)]
    X[event == "end",       (compartments) := list(+0L, +0L, +0L)]

    X[, event := NULL]
    X[, S := cumsum(S)]
    X[, I := cumsum(I)]
    X[, R := cumsum(R)]

    X <- X[time <= tmax]

    X
}

make_time_series_sis <- function(popn, params) {
    compartments <- c("S", "I")
    # tmax <- params$tmax
    tmax <- max(popn$Tdeath, na.rm = TRUE)

    popn2 <- popn[sdp == "progeny", .(Tinf, Tdeath)]
    N <- popn2[, .N]

    X <- rbind(
        data.table(event = factor("start", c("start", "infection", "removal", "end")),
                   time = 0.0),
        data.table(event = "infection", time = popn$Tinf[popn2$Tinf > 0]),
        data.table(event = "removal",   time = popn$Tdeath),
        data.table(event = "end",       time = tmax))

    X <- X[is.finite(time)]

    setorder(X, time)

    X[event == "start",     (compartments) := list(+N,  +0L)]
    X[event == "infection", (compartments) := list(-1L, +1L)]
    X[event == "removal",   (compartments) := list(+1L, -1L)]
    X[event == "end",       (compartments) := list(+0L, +0L)]

    X[, event := NULL]
    X[, S := cumsum(S)]
    X[, I := cumsum(I)]

    X <- X[time <= tmax]

    X
}
