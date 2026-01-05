source("censor_data.R")
source("discretise_time_sire.R")
source("utils.R")

# Tidy up ready for saving data to XML file
prepare_data <- function(popn, params) {
    message("Reformatting data for SIRE ...")

    {
        time_step        <- params$time_step
        compartments     <- params$compartments
        tmax             <- params$tmax
        traitnames       <- params$traitnames
        model_traits     <- params$model_traits
        trial_fe         <- params$trial_fe
        donor_fe         <- params$donor_fe
        txd_fe           <- params$txd_fe
        weight_fe        <- params$weight_fe
        weight1_fe       <- params$weight1_fe
        weight2_fe       <- params$weight2_fe
        use_weight       <- params$use_weight
        weight_is_nested <- params$weight_is_nested
        sim_new_data     <- params$sim_new_data
        censor           <- params$censor
        timings          <- params$timings
        pass_events      <- params$pass_events
    }
    
    data <- copy(popn)
    
    # Sire and dam need "."s for missing values
    data[, sire := fifelse(is.na(sire), ".", as.character(sire))]
    data[, dam  := fifelse(is.na(dam),  ".", as.character(dam))]
    
    
    
    # Fixed effects
    data[, trial_fe := fifelse(trial == 2, 1, 0, 0)]
    data[, donor_fe := fifelse(donor == 1, 1, 0, 0)]
    data[, txd_fe   := fifelse(trial == 2 & donor == 1, 1, 0, 0)]
    
    if ("weight" %in% names(data)) {
        rec_fn <- if (use_weight == "log") log_recentre else recentre
        if (weight_is_nested) {
            data[, `:=`(weight1_fe = 0, weight2_fe = 0)]
            data[trial == 1, weight1_fe := rec_fn(weight)]
            data[trial == 2, weight2_fe := rec_fn(weight)]
        } else {
            data[, weight_fe := 0]
            data[sdp == "progeny", weight_fe := rec_fn(weight)]
        }
    }
    
    # Initial compartment
    data[, initial_comp := compartments[donor + 1]]
    setcolorder(data, "initial_comp", after = "donor")
    
    # Tidy up BVs
    if (any(str_detect(names(data), "_g"))) {
        TBVs <- intersect(str_c(traitnames, "_g"), names(data))
        TBV1 <- str_c(str_sub(TBVs, 1, 1), "_g")
        setnames(data, TBVs, TBV1)
        # We don't really need more than 4 S.F.s
        data[, (TBV1) := map(.SD, signif, 4), .SDcols = TBV1]
    } else {
        TBVs <- NULL
        TBV1 <- NULL
    }
    
    
    # Timings ----
    
    # Set any timings not in pass_events to NA
    walk(setdiff(timings, pass_events), \(ev) {
        # Donors get to keep their Tinfs, but that's it
        data[, (ev) := fifelse(ev == "Tinf" & donor == 1, 0, NA)]
    })
    
    # If we've added a new column, group the timings together for ease of
    # checking the data file.
    setcolorder(data, timings, after = "initial_comp")
    
    
    # Censoring
    #
    # This needs to be done before discretising, but might need to be redone if
    # any of the timings are bunched together too much.
    if (sim_new_data != "no" && censor < 1) {
        walk(intersect(timings, names(data)), \(timing) {
            data[get(timing) >= tmax[str_c("t", trial)], (timing) := NA]
        })
    }
    
    
    # Discretise time
    
    # tidy up to 3 decimal places
    cols <- intersect(timings, names(data))
    data[, (cols) := map(.SD, round, 3), .SDcols = cols]
    
    if (time_step > 0) {
        data[, comp_status := discretise_time_sire(data, params)]
    } else {
        # Check that Tinf < Tinc < Tsym < Tdet, otherwise bump up each timing to
        # at least 0.5 more than the previous timing
        walk(seq_along(timings)[-1], \(i) {
            ti1 <- timings[[i - 1]]
            ti2 <- timings[[i]]
            data[get(ti1) >= get(ti2), (ti2) := get(ti1) + 0.5]
            
            # re-censor if anything has moved past tmax
            data[get(ti2) >= tmax[str_c("t", trial)], (ti2) := NA]
        })
    }
    
    
    # Convert columns to strings with "no" and "."
    fix_NAs <- function(x, alt = "no") fifelse(is.na(x), alt, x)
    
    cols <- c(timings, "group", TBV1)
    data[,                 (cols) := map(.SD, as.character),  .SDcols = cols]
    data[sdp != "progeny", (cols) := map(.SD, fix_NAs, "NA"), .SDcols = cols]
    data[sdp == "progeny", (cols) := map(.SD, fix_NAs, "no"), .SDcols = cols]
    
    # If we see "no" for an event, but a later event is not "no", it's actually
    # missing, so change to "."
    walk(rev(seq_along(timings)[-1]), \(i) {
        ti1 <- timings[[i - 1]]
        ti2 <- timings[[i]]
        data[get(ti1) == "no" & get(ti2) != "no" & Tdeath != tmax[str_c("t", trial)], (ti1) := "."]
    })
    
    # Donors inoculation event occurs outside of the observation times
    data[donor == 1, (timings[[1]]) := "no"]
    
    
    # Tidy up and remove unwanted columns (including trait PTs and EVs)
    cols_to_del <- intersect(
        c("GE", "status", "generation", "infected_by", "parasites", "donor", "weight",
          model_traits, str_c(traitnames, "_EV")),
        names(data))
    
    data[, (cols_to_del) := NULL]
    
    
    data
}
