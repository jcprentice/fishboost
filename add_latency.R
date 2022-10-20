# Add latency period (LP) to pop based on average from donors.

add_latency_donors <- function(pop) {
    # simulated data doesn't have a trial, so let it temporarily = 1
    has_trial <- "trial" %in% names(pop)
    if (!has_trial)
        pop[, trial := 1]

    LP <- pop[sdp == "progeny" & donor == 1,
              .(LP = mean(Tsym, na.rm = TRUE)),
              by = trial][, LP]
    message(" - Adding LP = ", round(LP[1], digits = 1), " days")
    pop[, estimated_Tinf := fifelse(donor == 1, 0, pmax(Tsym - LP[trial], 0))]

    # tidy up trial
    if (!has_trial)
        pop[, trial := NULL]
}


# Add latency period (LP) to pop based on average from donors, and extend to
# fill any gaps. Round times down.


add_latency_donors_plus <- function(pop) {
    # simulated data doesn't have a trial, so let it temporarily = 1
    has_trial <- "trial" %in% names(pop)
    if (!has_trial)
        pop[, trial := 1]

    LP_donors <- pop[sdp == "progeny" & donor == 1,
              .(LP = ceiling(mean(Tsym, na.rm = TRUE))),
              by = trial][, LP]
    message(" - Adding LP = ", LP_donors[1], " days")

    setkey(pop, group, Tsym, Trec)

    pop[!is.na(Trec) & !is.na(Tsym),
        LP := pmax(Tsym - cummax(shift(Trec - 1, n = 1, fill = 0)), 0),
        by = group]
    LP_extended <- pop[, ceiling(max(LP, na.rm = TRUE))]
    message(" - Extending LP by up to ", LP_extended, " days")

    pop[, LP := pmax(LP, LP_donors[trial])]
    pop[, estimated_Tinf := fifelse(donor == 1, 0, pmax(Tsym - LP, 0))]

    setkey(pop, id)

    # tidy up trial
    if (!has_trial)
        pop[, trial := NULL]
}


add_latency_biggest <- function(pop) {
    setkey(pop, group, Tsym, Trec)

    pop[sdp == "progeny" & !is.na(Trec) & !is.na(Tsym),
        LP := pmax(Tsym - cummax(shift(Trec - 1, n = 1, fill = 0)), 0),
        by = group]
    LP_extended <- pop[, ceiling(max(LP, na.rm = TRUE))]
    message(" - Extending LP by up to ", LP_extended, " days")

    pop[, LP := pmax(LP, LP_donors[trial])]
    pop[, estimated_Tinf := fifelse(donor == 1, 0, pmax(Tsym - LP, 0))]

    setkey(pop, id)

    # tidy up trial
    if (!has_trial)
        pop[, trial := NULL]
}


# Add latency period (LP) to pop. Calculate smallest LP necessary to cover all
# gaps in infection, apply to ALL individuals

add_latency_smallest <- function(pop) {
    pop1 <- pop[sdp == "progeny", .(group, Tsym, Trec)]
    setkey(pop1, group, Tsym, Trec)

    # Remove individuals who didn't get infected
    # tmp1 <- pop1[!is.na(Trec)]
    # tmp2 <- tmp1[, .(lp = pmax(Tsym - cummax(shift(Trec, n = 1, fill = 0)), 0)), by = group]
    # lp   <- tmp2[, max(lp)]

    lp <- pop1[!is.na(Trec)
    ][, .(lp = pmax(Tsym - cummax(shift(Trec, n = 1, fill = 0)), 0)), by = group
    ][, max(lp)]

    message(" - Adding LP = ", round(lp, digits = 1), " days")

    pop[, estimated_Tinf := pmax(Tsym - lp, 0)]
}



# Add latency period (LP) to pop. Calculate smallest LP necessary to cover all
# gaps in infection for each individual separately.

add_latency_individual <- function(pop) {
    pop1 <- pop[sdp == "progeny", .(id, group, Tsym, Trec)]
    setkey(pop1, group, Tsym, Trec)

    pop2 <- pop1[!is.na(Trec)][, .(id, Tinf2 = pmin(cummax(shift(Trec, n = 1, fill = 0)), Tsym)), by = group]
    pop3 <- merge(pop[, .(id)], pop2[, .(id, Tinf2)], by = "id", all = TRUE)

    pop[, estimated_Tinf := pop3$Tinf2]

    rng <- pop[, round(range(Tsym - Tinf, na.rm = TRUE), digits = 1)]

    message(" - Adding LPs in [", rng[1], ", ", rng[2], "] days")
}


