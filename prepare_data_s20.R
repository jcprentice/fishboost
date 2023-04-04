source("add_latency.R")


# Tidy up ready for saving data to XML file

prepare_data <- function(pop, params) {
    message("Tidying up trait and pop data...")

    model_type   <- params$model_type
    t_gap        <- params$t_gap
    pass_Tsym    <- params$pass_Tsym
    use_sire_21  <- params$use_sire_21
    compartments <- params$compartments
    traitnames   <- params$traitnames

    # quick fix as the SIR model doesn't have a separate Tsym
    if (model_type == "SIR") {
        pass_Tsym <- "Tinf"
    }

    # SIRE 2.0 wants Tinf, but what we give it is actually Tsym
    switch(pass_Tsym,
           "Tinf" = {
               message(" - passing 'Tinf'")
               timings <- c("Tinf", "Trec")
           }, "Tsym" = {
               message(" - passing 'Tsym'")
               timings <- c("Tsym", "Trec")
           }, "estimated_Tinf_from_donors" = {
               add_latency_donors(pop)
               timings <- c("estimated_Tinf", "Trec")
           }, "estimated_Tinf_from_donors_plus" = {
               add_latency_donors_plus(pop)
               timings <- c("estimated_Tinf", "Trec")
           }, "estimated_Tinf_per_individual" = {
               add_latency_individual(pop)
               timings <- c("estimated_Tinf", "Trec")
           }, stop(" - Didn't understand 'pass_Tsym = \"", pass_Tsym, "\"'")
    )

    # these are the columns we want for the data table
    sir <- c("susceptibility", "infectivity", "recoverability")
    used <- sir %in% traitnames
    traitnames_BV <- if (any(used)) {
        paste0(sir[sir %in% traitnames], "_BV")
    } else {
        NULL
    }
    cols <- c("id", "sdp", "group", "donor", timings, traitnames_BV)

    traits2 <- pop[sdp != "dam", ..cols]
    traits2[, id := paste0("Ind", seq.int(0L, .N - 1L))]

    # Whatever we got, call it "Tsym". Note: SIRE 2.0 calls it "It"
    setnames(traits2, timings[1], "Tsym")

    # tmax <- ceiling(max(traits2[, .(Tsym, Trec)], na.rm = TRUE))
    tmax <- ceiling(max(traits2$Trec, na.rm = TRUE))

    # handle conversion to discrete time
    if (t_gap > 0) {
        message(" - Applying Tgap = ", t_gap)
        tmax <- ceiling(tmax / t_gap) * t_gap
        N <- nrow(traits2)
        Tsym <- ceiling(traits2$Tsym / t_gap) * t_gap
        Trec <- ceiling(traits2$Trec / t_gap) * t_gap

        state <- rep(NA_character_, N)

        for (i in traits2[, .I[sdp == "progeny"]]) {
            # [S,0],[I,2],[I,20],[R,22]
            St <- 0; It <- Tsym[i]; Rt <- Trec[i]

            # handle S
            str <- character(0)
            # [S,0] only if neither It or Rt are 0
            if ((is.na(It) && is.na(Rt)) || (!is.na(It) && It > 0) || (is.na(It) && !is.na(Rt) && Rt > 0)) {
                str <- c(str, "[S,0]")
            }
            # [S, last val]
            if (!is.na(It) && It > St + t_gap) {
                str <- c(str, paste0("[S,", It - t_gap, "]"))
            } else if (is.na(It) && !is.na(Rt) && Rt > St + t_gap) {
                str <- c(str, paste0("[S,", Rt - t_gap, "]"))
            } else if (is.na(It) && is.na(Rt)) {
                str <- c(str, paste0("[S,", tmax, "]"))
            }
            # handle I
            if (!is.na(It)) {
                str <- c(str, paste0("[I,", It, "]"))
            }
            # [I, last val]
            if (!is.na(It) && !is.na(Rt) && Rt > It + t_gap) {
                str <- c(str, paste0("[I,", Rt - t_gap, "]"))
            } else if (!is.na(It) && is.na(Rt)) {
                str <- c(str, paste0("[I,", tmax, "]"))
            }
            # handle R
            if (!is.na(Rt)) {
                str <- c(str, paste0("[R,", Rt, "]"))
            }
            # [R, last val]
            if (!is.na(Rt) && Rt < tmax - t_gap) {
                str <- c(str, paste0("[R,", tmax, "]"))
            }

            state[i] <- paste(str, collapse = ",")
        }
        traits2[, "state"] <- state
    }

    # In SIRE 2.0 traits act inversely to SIRE 2.1
    if (use_sire_21 == FALSE) {
        traits2[, recoverability_BV := - recoverability_BV]
    }

    # Fix Tsym == Trec, set Tsym <- Tsym -1
    traits2[Tsym == Trec & Tsym > 1, Tsym := Tsym - 1]

    # Handle NA cases ----

    # SIRE 2.1 needs initial_comp, and also uses "NA" instead of "." for parent
    # values. Here donors have value 1, receptors have value 0, parents have no
    # value, so c("S", "E", "I", "R")[1, 0, NA] = ["E", "S", NA].
    traits2[, donor := compartments[donor + 1]]

    NA_char <- if (use_sire_21) "NA" else "."

    fix_NAs <- function(x, alt = "no") ifelse(is.na(x), alt, as.character(x))

    cols1 <- c(if (t_gap == 0) c("Tsym", "Trec") else "state",
               "group", traitnames_BV)

    traits2[,                 (cols1) := lapply(.SD, as.character),     .SDcols = cols1]
    traits2[sdp == "sire",    (cols1) := lapply(.SD, fix_NAs, NA_char), .SDcols = cols1]
    traits2[sdp == "progeny", (cols1) := lapply(.SD, fix_NAs, "no"),    .SDcols = cols1]

    if (t_gap > 0) {
        traits2[, c("Tsym", "Trec") := NULL]
        setcolorder(traits2, c(1:3, ncol(traits2)))
    }

    # Fix some problems for SIRE 2.1
    if (use_sire_21) {
        # Transition time for S->I for donors does not occur within the
        # observation period
        if (model_type %in% c("SIR", "SIR_res")) {
            traits2[donor == compartments[2], Tsym := "no"]
        } else {
            traits2[donor == compartments[2] & Tsym == "no", Tsym := "."]
        }

        # Shouldn't have any Tsym = 0
        traits2[Tsym == 0, Tsym := 1]

        # Fix missing Tsyms for secondary infectives
        traits2[donor == "S" & Tsym == "no" & Trec != "no", Tsym := "."]
    }

    # SIRE 2.0 doesn't want the column "donor"
    if (use_sire_21 == FALSE) {
        traits2[, donor := NULL]
    }

    # Remove any other unwanted columns
    traits2[, sdp := NULL]

    traits2
}
