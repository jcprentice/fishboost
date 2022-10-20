make_time_series <- function(pop, params) {
    switch(params$model_type,
           "SIR"  = make_time_series_sir(pop, params),
           "SEIR" = make_time_series_seir(pop, params),
           "SEIDR" = make_time_series_seidr(pop, params)
    )
}


make_time_series_seidr <- function(pop, params) {
    compartments <- c("S", "E", "I", "D", "R")
    # tmax <- params$tmax
    tmax <- max(pop$Trec, na.rm = TRUE)

    pop2 <- pop[sdp == "progeny", .(Tinf, Tinc, Tsym, Trec)]
    pop2[is.na(Tsym), Tsym := Trec]
    pop2[is.na(Tinc), Tinc := Tsym]
    pop2[is.na(Tinf), Tinf := Tinc]

    N <- pop2[, .N]

    X <- rbind(
        data.table(event = factor("start", c("start", "infection", "incubation", "detection", "recovery", "end")),
                   time = 0.0),
        data.table(event = "infection",  time = pop2$Tinf),
        data.table(event = "incubation", time = pop2$Tinc),
        data.table(event = "detection",  time = pop2$Tsym),
        data.table(event = "recovery",   time = pop2$Trec),
        data.table(event = "end", time = tmax))

    X <- X[is.finite(time)]

    setkey(X, time)
    #                                                S    E    I    D    R
    X[event == "start",      (compartments) := list(+N,  +0L, +0L, +0L, +0L)]
    X[event == "infection",  (compartments) := list(-1L, +1L, +0L, +0L, +0L)]
    X[event == "incubation", (compartments) := list(+0L, -1L, +1L, +0L, +0L)]
    X[event == "detection",  (compartments) := list(+0L, +0L, -1L, +1L, +0L)]
    X[event == "recovery",   (compartments) := list(+0L, +0L, +0L, -1L, +1L)]
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


make_time_series_sidr <- function(pop, params) {
    compartments <- c("S", "I", "D", "R")
    # tmax <- params$tmax
    tmax <- max(pop$Trec, na.rm = TRUE)

    pop2 <- pop[sdp == "progeny", .(Tinf, Tsym, Trec)]
    pop2[is.na(Tsym), Tsym := Trec]
    pop2[is.na(Tinf), Tinf := Tsym]

    N <- pop2[, .N]

    X <- rbind(
        data.table(event = factor("start", c("start", "infection", "detection", "recovery", "end")),
                   time = 0.0),
        data.table(event = "infection",  time = pop2$Tinf),
        data.table(event = "detection",  time = pop2$Tsym),
        data.table(event = "recovery",   time = pop2$Trec),
        data.table(event = "end", time = tmax))

    X <- X[is.finite(time)]

    setkey(X, time)
    #                                                S    I    D    R
    X[event == "start",      (compartments) := list(+N,  +0L, +0L, +0L)]
    X[event == "infection",  (compartments) := list(-1L, +1L, +0L, +0L)]
    X[event == "detection",  (compartments) := list(+0L, -1L, +1L, +0L)]
    X[event == "recovery",   (compartments) := list(+0L, +0L, -1L, +1L)]
    X[event == "end",        (compartments) := list(+0L, +0L, +0L, +0L)]

    X[, event := NULL]
    X[, S := cumsum(S)]
    X[, I := cumsum(I)]
    X[, D := cumsum(D)]
    X[, R := cumsum(R)]
    X[, ID := I + D]

    X <- X[time <= tmax]

    X
}


make_time_series_seir <- function(pop, params) {
    compartments <- c("S", "E", "I", "R")
    # tmax <- params$tmax
    tmax <- max(pop$Trec, na.rm = TRUE)

    pop2 <- pop[sdp == "progeny", .(Tinf, Tsym, Trec)]
    pop2[is.na(Tsym), Tsym := Trec]
    pop2[is.na(Tinf), Tinf := Tsym]

    N <- pop2[, .N]

    X <- rbind(
        data.table(event = factor("start", c("start", "infection", "incubation", "recovery", "end")),
                   time = 0.0),
        data.table(event = "infection",  time = pop2$Tinf),
        data.table(event = "incubation", time = pop2$Tsym),
        data.table(event = "recovery",   time = pop2$Trec),
        data.table(event = "end", time = tmax))

    X <- X[is.finite(time)]

    setkey(X, time)

    X[event == "start",      (compartments) := list(+N,  +0L, +0L, +0L)]
    X[event == "infection",  (compartments) := list(-1L, +1L, +0L, +0L)]
    X[event == "incubation", (compartments) := list(+0L, -1L, +1L, +0L)]
    X[event == "recovery",   (compartments) := list(+0L, +0L, -1L, +1L)]
    X[event == "end",        (compartments) := list(+0L, +0L, +0L, +0L)]

    X[, event := NULL]
    X[, S := cumsum(S)]
    X[, E := cumsum(E)]
    X[, I := cumsum(I)]
    X[, R := cumsum(R)]

    X <- X[time <= tmax]

    X
}


make_time_series_sir <- function(pop, params) {
    compartments <- c("S", "I", "R")
    # tmax <- params$tmax
    tmax <- max(pop$Trec, na.rm = TRUE)

    pop2 <- pop[sdp == "progeny", .(Tinf, Trec)]
    pop2[is.na(Tinf), Tinf := Trec]

    N <- pop2[, .N]

    X <- rbind(
        data.table(event = factor("start", c("start", "infection", "recovery", "end")),
                   time = 0.0),
        data.table(event = "infection",  time = pop2$Tinf),
        data.table(event = "recovery",   time = pop2$Trec),
        data.table(event = "end", time = tmax))

    X <- X[is.finite(time)]

    setkey(X, time)

    X[event == "start",     (compartments) := list(+N,  +0L, +0L)]
    X[event == "infection", (compartments) := list(-1L, +1L, +0L)]
    X[event == "recovery",  (compartments) := list(+0L, -1L, +1L)]
    X[event == "end",       (compartments) := list(+0L, +0L, +0L)]

    X[, event := NULL]
    X[, S := cumsum(S)]
    X[, I := cumsum(I)]
    X[, R := cumsum(R)]

    X <- X[time <= tmax]

    X
}

make_time_series_sis <- function(pop, params) {
    compartments <- c("S", "I")
    # tmax <- params$tmax
    tmax <- max(pop$Trec, na.rm = TRUE)

    pop2 <- pop[sdp == "progeny", .(Tinf, Trec)]
    N <- pop2[, .N]

    X <- rbind(
        data.table(event = factor("start", c("start", "infection", "recovery", "end")),
                   time = 0.0),
        data.table(event = "infection",  time = pop$Tinf),
        data.table(event = "recovery",   time = pop$Trec),
        data.table(event = "end", time = tmax))

    X <- X[is.finite(time)]

    setkey(X, time)

    X[event == "start",     (compartments) := list(+N,  +0L)]
    X[event == "infection", (compartments) := list(-1L, +1L)]
    X[event == "recovery",  (compartments) := list(+1L, -1L)]
    X[event == "end",       (compartments) := list(+0L, +0L)]

    X[, event := NULL]
    X[, S := cumsum(S)]
    X[, I := cumsum(I)]

    X <- X[time <= tmax]

    X
}
