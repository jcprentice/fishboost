apply_time_step <- function(data, params) {
    {
        model_type    <- params$model_type
        tmax          <- params$tmax
        time_step     <- params$time_step
        ntraits       <- params$ntraits
        compartments  <- params$compartments
        use_parasites <- params$use_parasites
    }
    
    # Fishboost uses time_step = 1, I suppose I might want to overrule it though
    if (params$use_fb_data) {
        time_step <- 1
    }
    
    # This function finds the smallest time step >= x, ensuring x is numeric
    ceil_ts <- function(x, t = 1) ceiling(as.numeric(x) / t) * t
    
    # this is [0, ..., 104 or 160]
    all_times <- list()
    for (i in seq_along(tmax)) {
        # need to round up to nearest time_step to include all values
        all_times[[i]] <- seq(0, ceil_ts(tmax[i], time_step), by = time_step)
    }
    
    # If we only have 1 trial we run into problems indexing trial 2...
    ntrials <- data[sdp == "progeny", length(unique(trial))]
    
    # TODO: test if this works properly for non-SEIDR models
    unknown_states <- switch(
        model_type,
        "SEIDR" = c("S|E|I|D", "E|I|D", # missing Tsym
                    "S|E|I", "E|I",     # SEI
                    "S|E", "E",         # SE
                    "S", "S",           # S
                    "E|I", "E|I",       # EI
                    "I", "I"),          # I
        "SIDR" = c("S|I|D", "I|D", # missing Tsym
                   "S|I", "I",     # Tsym or SEI
                   "S", "I",       # Tsym or SE
                   "S", "S"),      # Tsym or S
        "SEIR" = c("S|E|I", "E|I",  # missing Tsym
                   "S|E|I", "E|I",  # 
                   "S|E|I", "E|I",  # 
                   "S|E|I", "E|I"), # 
        "SIR" = c("S|I", "I", # missing Tsym
                  "S|I", "I",
                  "S|I", "I",
                  "S|I", "I")
    )
    us = matrix(unknown_states, ncol = 2, byrow = TRUE,
                dimnames = list(c("missing", "SEI", "SE", "S", "EI", "I"),
                                c("R", "D")))
    
    # This function clips everything but the head and tail, and the union
    # prevents duplication
    head_tail <- function(x) {
        a <- x[1]; b <- x[length(x)]
        if (a == b) a else c(a, b)
    }
    
    comp_status <- character(nrow(data))
    
    # for (i in data[group == 2, id]) {
    for (i in seq_len(nrow(data))) {
        if (data$sdp[i] != "progeny") {
            comp_status[i] <- "NA"
            next
        }
        
        Tsym <- data$Tsym[i]
        Trec <- data$Trec[i]
        trial <- if (ntrials == 1) 1 else data$trial[i]
        times <- all_times[[trial]]
        tr_max <- times[length(times)]
        init_comp <- data$initial_comp[i]
        is_donor <- init_comp == compartments[2]
        has_parasites <- data$parasites[i]
        
        # if (Trec %in% c("no", ".") && !parasites) {
        #     message(glue("found issue at i={i}"))
        #     break
        # }
        # 
        # message(glue("i={i}: Tsym={Tsym}, Trec={Trec}, donor={is_donor}, parasites={has_parasites}"))
        
        if (Tsym == "no") {
            if (use_parasites == "") {
                # Not using parasites column, just give the details we know
                status <- c(paste0("[", init_comp, ":0]"),
                            paste0("[", us["missing", 1 + is_donor], ":", head_tail(times[-1]), "]"))
            } else if (has_parasites) {
                # Parasites present, ensure fish ends up infected
                status <- c(paste0("[", init_comp, ":0]"),
                            paste0("[", us["SEI", "R"], ":", times[c(2, length(times) - 1)], "]"),
                            paste0("[", us["EI", 1], ":", tr_max, "]"))
            } else if (is_donor) {
                # This is the one donor without parasites that we force to start Exposed
                status <- c(paste0("[", init_comp, ":0]"), # E at start
                            paste0("[", us["EI", "D"], ":", times[c(2, length(times) - 1)], "]"), # E|I throughout
                            paste0("[", us["I", "D"], ":", tr_max, "]")) # must be infected by the end
            } else {
                # No parasites, downgrade possible infection status
                status <- c(paste0("[S:0]"), # S at start
                            paste0("[", us[use_parasites, "R"], ":", c(times[2], tr_max), "]")) # S|E|I throughout
            }
        } else if (Tsym == ".") {
            # We didn't see any symptoms, but it could still be infected
            v <- unknown_states[1 + is_donor] # "S|E|I|D", "E|I|D"
            t1 <- head_tail(which(times < as.numeric(Trec))[-1])
            t2 <- head_tail(which(times >= as.numeric(Trec)))
            
            status <- c(paste0("[", init_comp, ":0]"),
                        if (length(t1) > 0) paste0("[", v, ":", times[t1], "]") else NULL,
                        if (length(t2) > 0) paste0("[R:", times[t2], "]") else NULL)
        } else {
            # Tsym has a value, so may need to use numeric versions ceil'd to
            # nearest time step
            Tsym1 <- ceil_ts(Tsym, time_step)
            Trec1 <- if (Trec == "no") NA_real_ else ceil_ts(Trec, time_step)
            
            v <- unknown_states[3 + is_donor] # "S|E|I", "E|I"
            # Need to handle odd case where Tsym = 1
            ti1 <- which(times < Tsym1)
            t1 <- if (length(ti1) < 2) NULL else head_tail(ti1[-1])
            if (Trec == "no") {
                t2 <- head_tail(which(times >= Tsym1))
                t3 <- t4 <- NULL
            } else if (Tsym1 == Trec1) {
                t2 <- NULL
                t3 <- head_tail(which(times >= Trec1))
                t4 <- Tsym1
            } else {
                t2 <- head_tail(which(times >= Tsym1 & times < Trec1))
                t3 <- head_tail(which(times >= Trec1))
                t4 <- NULL
            }
            
            status <- c(paste0("[", init_comp, ":0]"),
                        if (length(t1) > 0) paste0("[", v, ":", times[t1], "]") else NULL,
                        if (length(t2) > 0) paste0("[D:", times[t2], "]") else NULL,
                        if (length(t4) > 0) paste0("[D|R:", Tsym, "]") else NULL,
                        if (length(t3) > 0) paste0("[R:", times[t3], "]") else NULL)
        }
        
        # Collapse the matrix and remove multiples of ,,,,,
        status_c <- paste0(status, collapse = ",")
        comp_status[i] <- gsub(",+", ",", status_c)
    }
    
    comp_status
}
