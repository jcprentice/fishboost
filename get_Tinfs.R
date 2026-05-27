{
    library(stringr)
    library(data.table)
    library(purrr)
}

get_Tinfs <- function(dataset = "fb-test", scen = 1, rep = 1, nsamples = 0) {
    if (FALSE) {
        dataset <- "fb-test"; scen <- 1; rep <- 1; nsamples <- 0;
    }

    f <- str_glue("datasets/{dataset}/data/scen-{scen}-{rep}-out/summary_inf.rds")
    popn <- readRDS(f)$popn

    states <- unique(popn$state)
    if (nsamples > 0) {
        states <- c(rep(states, nsamples %/% length(states)),
                    sample(states, nsamples %% length(states))) |>
            sort()
    }

    popn[sdp == "progeny" & state %in% states, .(state, id, Tinf)] |>
        dcast(id ~ state, value.var = "Tinf")
}

get_median_Tinfs <- function(params) {
    if (FALSE) {
        dataset <- "fb-test"; scen <- 1; rep <- 1
    } else {
        dataset <- params$dataset
        scen <- params$scen
        rep <- params$replicate
    }

    f <- str_glue("datasets/{dataset}/data/scen-{scen}-{rep}-out/summary_inf.rds")
    popn <- readRDS(f)$popn

    popn[is.na(Tinf), Tinf := Inf]
    popn[sdp == "progeny", .(Tinf = median(Tinf)), id][, Tinf]
}
