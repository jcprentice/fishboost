censor_data <- function(popn, params) {
    tmax <- params$tmax

    popn2 <- copy(popn)

    popn2[Tsign  %notin% c("NA", "no"), Tsign  := fifelse(as.numeric(Tsign)  >= tmax[str_c("t", trial)], "no", Tsign)]
    popn2[Tdeath %notin% c("NA", "no"), Tdeath := fifelse(as.numeric(Tdeath) >= tmax[str_c("t", trial)], "no", Tdeath)]

    popn2
}
