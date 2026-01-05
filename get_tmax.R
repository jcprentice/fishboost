# Find tmax for each trial
get_tmax <- function(popn, params) {
    censor <- params$censor
    time_step <- params$time_step
    sim_new_data <- params$sim_new_data
    
    if (sim_new_data == "no") {
        return(c(t1 = 104, t2 = 160))
    }
    
    x1 <- popn[sdp == "progeny", .(id, trial, Tdeath)]
    
    x1[, Tdeath := fifelse(is.na(Tdeath),
                           max(Tdeath, na.rm = TRUE),
                           Tdeath),
       trial]
    
    x2 <- x1[, quantile(Tdeath, censor, na.rm = TRUE), trial]
    tmax <- setNames(x2$V1, str_c("t", x2$trial))
    
    # Want the nearest time step if time is measured in discrete intervals
    if (time_step > 0) {
        tmax <- ceiling(tmax / time_step) * time_step
    }
    
    tmax
}
