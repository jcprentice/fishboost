# Find tmax for each trial
get_tmax <- function(pop, params) {
    censor <- params$censor
    time_step <- params$time_step
    
    tmax <- if (TRUE) {
        # Find tmax for based on all individuals, including uninfected
        (pop[sdp == "progeny", .(trial, Trec)]
         # [is.na(Trec), Trec := Inf]
         [, quantile(Trec, censor, na.rm = TRUE), trial]
         [, V1])
    } else {
        # Find tmax based only on infected individuals
        pop[sdp == "progeny", quantile(Trec, censor, na.rm = TRUE), trial][, V1]
    }
    
    # Want the nearest time step if time is measured in discrete intervals
    if (time_step > 0) {
        tmax <- ceiling(tmax / time_step) * time_step
    }
    
    tmax
}