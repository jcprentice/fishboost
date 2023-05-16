library(data.table)
library(glue)

calc_fixed_effects <- function(dataset = "fb-parasites4",
                               scenario = 3) {
    
    files <- load(glue("results/{dataset}/scen-{scenario}-1.RData"))
    
    x <- parameter_estimates[startsWith(parameter, "trial") |
                                 startsWith(parameter, "donor") |
                                 startsWith(parameter, "txd"),
                             .(parameter, mean)]
    FEs <- x$mean
    names(FEs) <- x$parameter
    rm(x)
    
    # We need this to ensure donors have been properly allocated
    data <- fread(glue("data/{dataset}/scen-{scenario}-1-data.tsv"))
    data <- data[sdp == "progeny", .(group, trial_fe, donor_fe, txd_fe)]
    
    ptrial <- data[, mean(trial_fe)]
    pdonor <- data[, mean(donor_fe)]
    ptxd <- data[, mean(txd_fe)]
    rm(data)
    
    
    vals <- data.table(donor = rep(c("donor", "recipient"), 2),
                       trial = rep(1:2, each = 2),
                       lat = numeric(4), inf = numeric(4), det = numeric(4), rec = numeric(4))
    
    for (i in 1:4) for (trait in c("lat", "inf", "det", "rec")) {
        vals[i, (trait) := ((donor == "donor") - pdonor) * FEs[[paste0("donor_", substr(trait, 1 ,1))]] + 
                 ((trial == 2) - ptrial) * FEs[[paste0("trial_", substr(trait, 1 ,1))]] +
                 ((trial == 2 & donor == "donor") - ptxd) * FEs[[paste0("txd_", substr(trait, 1 ,1))]]]
    }
    vals
}

vals <- calc_fixed_effects()