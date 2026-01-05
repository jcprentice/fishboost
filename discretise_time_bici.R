discretise_time_bici <- function(popn, params) {
    {
        sim_new_data  <- params$sim_new_data
        model_type    <- params$model_type
        tmax          <- params$tmax
        time_step     <- params$time_step
        compartments  <- params$compartments
        fix_fb_data   <- params$fix_fb_data
        pass_events   <- params$pass_events
        timings       <- params$timings
    }
    
    popn2 <- copy(popn)
    
    # Fishboost uses time_step = 1, I suppose I might want to overrule it though
    if (sim_new_data == "no") {
        time_step <- 1
    }
    
    # This function finds the smallest time step >= x, ensuring x is numeric
    ceil_ts <- function(x) ceiling(x / time_step) * time_step
    
    # this is [0, ..., 104 or 160]
    all_times <- map(tmax, ~ seq(0, ceil_ts(.x), by = time_step))
    
    # This function clips everything but the head and tail, and the union
    # prevents duplication
    head_tail <- function(x) {
        a <- x[[1L]]
        b <- x[[length(x)]]
        if (a == b) a else c(a, b)
    }
    
    n_events <- length(timings)
    
    clip <- function(x, a, b) x[a <= x & x <= b]
    
    # Generate appropriate string (e.g. "S|E|I") for a given t
    possible_states <- function(t, event_times) {
        fwd <- nafill(c(0, event_times), "locf")
        bwd <- nafill(c(event_times, Inf), "nocb")
        str_flatten(compartments[fwd <= t & t < bwd], "|") # FIXME: t <= bwd or t < bwd?
    }
    
    final_comp <- compartments[[length(compartments)]]
    
    # We haven't added unused events
    popn2[, (setdiff(timings, pass_events)) := NA_real_]
    
    map(seq_len(nrow(popn2)), \(i) {
        # i <- 60
        if (popn2$sdp[i] != "progeny") return(NULL)
        
        # popn2[id == 1, .(id, Tinf, Tinc, Tsym, Tdeath, parasites, initial_comp, comp_status)]
        
        tr_max <- tmax[[str_c("t", popn2$trial[[i]])]]
        init_comp <- compartments[[popn2$donor[[i]] + 1L]]
        
        if (init_comp == final_comp) {
            return(data.frame(id = i, t = c(0, tr_max), DS = final_comp))
        }
        
        is_donor <- popn2$donor[[i]] != 0
        has_parasites <- "parasites" %in% fix_fb_data && popn2$parasites[[i]]
        event_times <- popn2[i, ..timings] |> unlist() |> ceil_ts()
        
        # No Tdeath means no Tdeath, not missing
        if (is.na(event_times[["Tdeath"]])) {
            event_times[["Tdeath"]] <- tr_max + time_step
            # If no Tdeath *and* no Tsym, check for parasites
            if (is.na(event_times[["Tsym"]] && has_parasites == FALSE)) {
                event_times[["Tsym"]] <- tr_max + time_step
            }
        }
        
        # mget(c("trial", "tr_max", "init_comp", "is_donor", "has_parasites", "event_times"))
        
        times1 <- c(0, time_step, event_times - time_step, event_times,
                    event_times + time_step, tr_max - time_step, tr_max) |>
            sort() |> unique() |> clip(0, tr_max)
        
        
        states <- map_chr(times1, possible_states, event_times)
        
        # Fix start and end
        states[[1]] <- init_comp
        ns <- length(states)
        if (is_donor) {
            states <- str_remove(states, "S\\|")
        } else if (has_parasites) {
            states[[ns]] <- str_remove(states[[ns]], "S\\|")
        }
        
        # Remove runs > 3
        idxs <- map_lgl(seq_along(states), \(idx) {
            if (idx == 1 || idx == ns) TRUE
            else states[[idx - 1]] != states[[idx + 1]]
        })
        
        data.frame(ID = i, t = times1[idxs], DS = states[idxs])
    }) |> rbindlist()
}
