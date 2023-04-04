censor_data <- function(data, params) {
    tmax <- params$tmax

    data[!Tsym %in% c("NA", "no"), Tsym := fifelse(as.numeric(Tsym) >= tmax[trial], "no", Tsym)]
    data[!Trec %in% c("NA", "no"), Trec := fifelse(as.numeric(Trec) >= tmax[trial], "no", Trec)]
}
