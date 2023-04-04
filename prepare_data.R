source("censor_data.R")
source("apply_time_step.R")

# Tidy up ready for saving data to XML file
prepare_data <- function(pop, params) {
    message("Reformatting data for SIRE...")

    {
        model_type    <- params$model_type
        timings       <- params$timings
        time_step     <- params$time_step
        pass_Tsym     <- params$pass_Tsym
        compartments  <- params$compartments
        tmax          <- params$tmax
        traitnames    <- params$traitnames
        trial_fe      <- params$trial_fe
        donor_fe      <- params$donor_fe
        txd_fe        <- params$txd_fe
        use_fb_data   <- params$use_fb_data
        censor        <- params$censor
        pass_events   <- params$pass_events
        use_parasites <- params$use_parasites
    }
    
    # pass the last N >= 2 events, but force N=2 if using FB data
    pass_events <- if (use_fb_data) 2 else max(pass_events, 2)
    
    timings <- tail(timings, pass_events)
    
    # these are the columns we want for the data table
    sir <- c("susceptibility", "infectivity", "recoverability")
    sir1 <- c("s", "i", "r")
    used <- sir %in% traitnames
    if (!use_fb_data && any(used)) {
        traitnames_BV <- paste0(sir[used], "_BV")
        trait_cols <- paste0(sir1[used], "_a")
    } else {
        traitnames_BV <- trait_cols <- NULL
    }

    # We need "trial" in data, so create it if it's not already there
    if (!"trial" %in% names(pop)) {
        pop[, trial := 0L]
    }
    
    # This corresponds to the order we're expecting to write them in
    cols <- c("id", "sire", "dam", "sdp", "trial", "group", "donor", "parasites",
              timings, traitnames_BV)
    
    # Make a copy of all the data
    data <- pop[, ..cols] # pop[sdp != "dam", ..cols] # exclude dams
    
    # data[, id := 1:.N] # paste0("Ind", seq.int(0L, .N - 1L))]

    # Whatever we got, call it "Tsym". Note: SIRE 2.0 calls it "It"
    setnames(data, timings[length(timings) - 1], "Tsym")

    # Traitnames go by s_a, i_a, r_a
    if (!use_fb_data && any(used)) {
        setnames(data, traitnames_BV, trait_cols)
    }
    
    
    # Fixed effects ----

    # 0 = {parents, trial 1}, 1 = trial 2
    if (trial_fe != "") {
        data[, trial_fe := fifelse(trial == 2, 1, 0, 0)]
    }
    
    # 0 = {parents, recipients}, 1 = donors
    if (donor_fe != "") {
        data[, donor_fe := fifelse(donor == 1, 1, 0, 0)]
    }

    # 0 = {parents, trial 1, trial 2 recipients}, 1 = trial 2 donors
    if (txd_fe != "") {
        data[, txd_fe := fifelse(trial == 2 & donor == 1, 1, 0, 0)]
    }
    

    # Handle time collisions ----
    
    if (time_step == 0) {
        # With continuous time, adjacent timings can't be equal, else dt = 0.
        
        # data[Tsym == Trec & Tsym < 1,  Trec := Trec + 0.5]
        # data[Tsym == Trec & Tsym >= 1, Tsym := Tsym - 0.5]
        data[Tsym == 0, `:=`(Tsym = 0.5, Trec = max(Trec, 0.5))]
        
        for (i in seq_along(timings)[-1]) {
            # get(timings[i]) is basically "Tsym" etc. If they're too close,
            # nudge the next one along up by 0.5
            data[get(timings[i-1]) >= get(timings[i]) - 0.5,
                 (timings[i]) := get(timings[i-1]) + 0.5]
        }
        # Just in case tmax has increased from this nudging
        # tmax <- data[sdp == "progeny", max(c(Tsym, Trec), na.rm = TRUE), trial][, V1]
        tmax <- get_tmax(data, params)
        params$tmax <- tmax
    } else {
        # This is already taken care of by SIRE for discrete time
    }

    # Don't really want or need higher accuracy than this
    for (t in timings) {
        data[, (t) := round(get(t), 3)]
    }
    for (t in trait_cols) {
        data[, (t) := signif(get(t), 4)]
    }
    
    # Handle NA cases ----
    
    # SIRE 2.1 needs initial_comp, and also uses "NA" instead of "." for parent
    # values. Here donors have value 1, receptors have value 0, parents have no
    # value, so c("S", "E", "I", "R")[0, 1, NA] => ["S", "E", NA].
    data[, initial_comp := compartments[donor + 1]]

    fix_NAs <- function(x, alt = "no") ifelse(is.na(x), alt, as.character(x))

    cols1 <- c(timings, "group", trait_cols)

    data[,                 (cols1) := lapply(.SD, as.character),  .SDcols = cols1]
    data[sdp != "progeny", (cols1) := lapply(.SD, fix_NAs, "NA"), .SDcols = cols1]
    data[sdp == "progeny", (cols1) := lapply(.SD, fix_NAs, "no"), .SDcols = cols1]

    cols2 <- c("sire", "dam")
    data[,                 (cols2) := lapply(.SD, as.character), .SDcols = cols2]
    data[sdp != "progeny", (cols2) := lapply(.SD, fix_NAs, "."), .SDcols = cols2]
    
    # In SIRE 2.1, transition time for S->E/I for donors does not occur within
    # the observation period, need to be careful with SIR models
    if (model_type %in% c("SIR", "SIR_res")) {
        data[initial_comp == compartments[2], Tsym := "no"]
    } else if ("Tinf" %in% timings) {
        data[initial_comp == compartments[2], Tinf := "no"]
    }
    
    setcolorder(data, c(1:6, ncol(data)))
    
    # Censoring ----
    
    if (use_fb_data) {
        # Correct for cases when Trec < end (i.e. the individual died from
        # infection) but Tsym is NA. Here Tsym is missing data, so needs to be
        # replaced with "."
        data[Tsym == "no" & !Trec %in% c("NA", "no") & Trec != tmax[trial], Tsym := "."]
        
        # No symptoms but parasites=T means it *could* have been infected. Since
        # the specificity = 1 and the sensitivity < 1, parasites=F means it was
        # definitely not infected).
        if (use_parasites != "") {
            # any donor (except 908 / 1830) with no symptoms and no parasites
            id_to_keep <- data[group == 71, id[1]]
            data[initial_comp == "E" & Tsym == "no" & !parasites & id != id_to_keep,
                 `:=`(initial_comp = "S", donor_fe = 0, txd_fe = 0)]
        }
        
        # In FB, Trecs at the trial ends are censored, not real Trecs, so set
        # this to "no"
        # data[Trec == as.character(tmax[trial]), Trec := "no"]
    } else if (censor < 1) {
        data <- censor_data(data, params)
    }

    # Discrete time ----
    if (time_step > 0) {
        data[, comp_status := apply_time_step(data, params)]
        # data[, intersect(timings, names(data)) := NULL]
    }
    
    # Tidy up ----
    data[, c("donor", "parasites") := NULL]

    data
}
