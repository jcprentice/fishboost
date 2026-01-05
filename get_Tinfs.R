{
    library(stringr)
    library(data.table)
    library(purrr)
}

get_Tinfs <- function(dataset = "fb-test", scen = 1, rep = 1, nsamples = 0) {
    # dataset <- "fb-test"; scen <- 1; rep <- 1; nsamples <- 0;
    
    f <- str_glue("datasets/{dataset}/data/scen-{scen}-{rep}-out/etc_inf.rds")
    popn <- readRDS(f)$popn
    
    states <- unique(popn$state)
    if (nsamples > 0) {
        states <- c(rep(states, nsamples %/% length(states)),
                    sample(states, nsamples %% length(states))) |>
            sort()
    }
    
    popn[state %in% states, .(state, id, Tinf)] |>
        dcast(id ~ state, value.var = "Tinf")
}

get_median_Tinfs <- function(dataset = "fb-test", scen = 1, rep = 1) {
    # dataset <- "fb-test"; scen <- 1; rep <- 1
    
    f <- str_glue("datasets/{dataset}/data/scen-{scen}-{rep}-out/etc_inf.rds")
    popn <- readRDS(f)$popn
    
    popn[is.na(Tinf), Tinf := Inf]
    popn[, .(Tinf = median(Tinf)), id]
}
