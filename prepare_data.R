source("apply_time_step.R")

# Tidy up ready for saving data to XML file
prepare_data <- function(pop, params) {
    message("Tidying up trait and pop data...")

    {
        model_type      <- params$model_type
        timings         <- params$timings
        time_step       <- params$time_step
        clear_Tsym      <- params$clear_Tsym
        pass_Tsym       <- params$pass_Tsym
        use_sire_21     <- params$use_sire_21
        compartments    <- params$compartments
        traitnames      <- params$traitnames
        trial_fe        <- params$trial_fe
        donor_fe        <- params$donor_fe
        txd_fe          <- params$txd_fe
        use_fb_data     <- params$use_fb_data
        ignore_donors   <- params$ignore_donors
        pass_events     <- params$pass_events
    }

    # pass the last N events, but force N=2 if using FB data
    pass_events <- max(pass_events, 2)
    if (use_fb_data) {
        pass_events <- 2
    }
    timings <- tail(timings, pass_events)
    
    # these are the columns we want for the data table
    sir <- c("susceptibility", "infectivity", "recoverability")
    used <- sir %in% traitnames
    traitnames_BV <- if (any(used) && !use_fb_data) {
        paste0(sir[used], "_BV")
    } else {
        c()
    }
    
    # easiest just to create trial and txd if they doesn't exist
    if (!"trial" %in% names(pop)) {
        pop[, trial := 0L]
    }
    if (!"txd" %in% names(pop)) {
        pop[, txd := 0L]
    }
    
    cols <- c("id", "sdp", "group", "trial", "donor", "txd",
              timings, traitnames_BV)


    # dams aren't included
    data <- pop[sdp != "dam", ..cols]
    data[, id := paste0("Ind", seq.int(0L, .N - 1L))]

    # Whatever we got, call it "Tsym". Note: SIRE 2.0 calls it "It"
    setnames(data, timings[length(timings) - 1], "Tsym")

    if (clear_Tsym) {
        data[, Tsym := NA]
    }

    # tmax <- ceiling(max(data[, .(Tsym, Trec)], na.rm = TRUE))
    tmax <- ceiling(max(data$Trec, na.rm = TRUE))

    
    # Fix Tsym == Trec, set Tsym <- Tsym -1
    data[Tsym == Trec & Tsym == 1, Trec := Trec + 0.5]
    data[Tsym == Trec & Tsym > 1,  Tsym := Tsym - 0.5]


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


    
    
    # Fix unreliable donors
    # data[donor == 1 & Tsym > ignore_donors, Tsym := NA]
    
    

    # Handle NA cases ----

    # SIRE 2.1 needs initial_comp, and also uses "NA" instead of "." for parent
    # values. Here donors have value 1, receptors have value 0, parents have no
    # value, so c("S", "E", "I", "R")[0, 1, NA] => ["S", "E", NA].
    data[, donor := compartments[donor + 1]]

    fix_NAs <- function(x, alt = "no") ifelse(is.na(x), alt, as.character(x))

    cols1 <- c(timings, "group", traitnames_BV)

    data[,                 (cols1) := lapply(.SD, as.character),  .SDcols = cols1]
    data[sdp == "sire",    (cols1) := lapply(.SD, fix_NAs, "NA"), .SDcols = cols1]
    data[sdp == "progeny", (cols1) := lapply(.SD, fix_NAs, "no"), .SDcols = cols1]

    
    # In SIRE 2.1, transition time for S->E/I for donors does not occur within
    # the observation period, need to be careful with SIR models
    if (model_type %in% c("SIR", "SIR_res")) {
        data[donor == compartments[2], Tsym := "no"]
    } else {
        if ("Tinf" %in% timings) {
            data[donor == compartments[2], Tinf := "no"]
        }
    }
    
    # Fix certain cases with the FB data:
    if (use_fb_data) {
        trial_ends <- c(104, 160)
        
        # Correct for cases when Trec < end (i.e. the individual died from
        # infection) but Tsym is NA. Here there should be a Tsym, but it was
        # unobserved
        data[Tsym == "no" & !Trec %in% c("NA", "no") & as.numeric(Trec) < trial_ends[trial], Tsym := "."]
        
        # In FB, Trecs at the trial ends are not real Trecs, so set this to "no"
        data[Trec == as.character(trial_ends[trial]), Trec := "no"]
    }

    if (time_step > 0) {
        data[, "Tsym"] <- apply_time_step(data, tmax, params)
        setnames(data, "Tsym", "comp_status")
        data[, Trec := NULL]
    } else {
        # Shouldn't have any Tsym = 0
        data[Tsym == 0, Tsym := 1]
        
        # Fix missing Tsyms for secondary infectives
        data[donor == "S" & Tsym == "no" & Trec != "no", Tsym := "."]
    }

    # Remove any other unwanted columns (note: SIRE 2.0 doesn't want "donor")
    data[, c("sdp", "trial", "txd") := NULL]

    data
}
