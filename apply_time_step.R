apply_time_step <- function(data, tmax, params) {
    {
        model_type   <- params$model_type
        ntotal       <- params$ntotal
        nsires       <- params$nsires
        nprogeny     <- params$nprogeny
        time_step    <- params$time_step
        ntraits      <- params$ntraits
        compartments <- params$compartments
    }
    
    # Fishboost uses time_step = 1, I suppose I might want to overrule it though
    if (params$use_fb_data) {
        time_step <- 1
    }
    
    # this is [0, ..., 104 or 160]
    times <- seq(0, tmax, by = time_step)
    
    # TODO: test if this works properly for non-SEIDR models
    unknown_states <- switch(
        model_type,
        "SEIDR" = c("S|E|I|D", "E|I|D", "S|E|I", "E|I"),
        "SIDR" = c("S|I|D", "I|D", "S|I", "I"),
        "SEIR" = c("S|E|I", "E|I", "S|E|I", "E|I"),
        "SIR" = c("S|I", "I", "S|I", "I")
    )
    
    # This function clips everything but the head and tail, and the union
    # prevents duplication
    head_tail <- function(x) union(head(x, 1), tail(x, 1))
    
    comp_status <- character(nsires + nprogeny)
    
    for (i in seq_len(nrow(data))) {
        if (data[i, sdp != "progeny"]) {
            comp_status[i] <- "NA"
            next()
        }
        
        Tsym <- data[i, Tsym]
        Trec <- data[i, Trec]
        donor <- data[i, donor]
        is_donor <- donor == compartments[2]
        
        if (Tsym == "no") {
            t1 <- head_tail(seq_along(times))
            status <- paste0("[", donor, ":", times[t1], "]")
        } else if (Tsym == ".") {
            v <- unknown_states[1 + is_donor]
            t1 <- head_tail(which(times < as.numeric(Trec))[-1])
            t2 <- head_tail(which(times >= as.numeric(Trec)))
            
            status <- c(glue("[{donor}:0]"),
                        if (length(t1) > 0) paste0("[", v, ":", times[t1], "]") else NULL,
                        if (length(t2) > 0) paste0("[R:", times[t2], "]") else NULL)
        } else {
            v <- unknown_states[3 + is_donor]
            t1 <- head_tail(which(times < as.numeric(Tsym))[-1])
            if (Trec == "no") {
                t2 <- head_tail(which(times >= as.numeric(Tsym)))
                t3 <- c()
            } else {
                t2 <- head_tail(which(times >= as.numeric(Tsym) & times < as.numeric(Trec)))
                t3 <- head_tail(which(times >= as.numeric(Trec)))
            }
            status <- c(glue("[{donor}:0]"),
                        if (length(t1) > 0) paste0("[", v, ":", times[t1], "]") else NULL,
                        if (length(t2) > 0) paste0("[D:", times[t2], "]") else NULL,
                        if (length(t3) > 0) paste0("[R:", times[t3], "]") else NULL)
        }
        
        # Collapse the matrix and remove multiples of ,,,,,
        status_c <- paste0(status, collapse = ",")
        comp_status[i] <- gsub(",+", ",", status_c)
    }
    
    comp_status
}

