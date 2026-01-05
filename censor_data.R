censor_data <- function(popn, params) {
    tmax <- params$tmax

    popn2 <- copy(popn)

    popn2[Tsym   %notin% c("NA", "no"), Tsym   := fifelse(as.numeric(Tsym)   >= tmax[str_c("t", trial)], "no", Tsym)]
    popn2[Tdeath %notin% c("NA", "no"), Tdeath := fifelse(as.numeric(Tdeath) >= tmax[str_c("t", trial)], "no", Tdeath)]

    popn2
}
