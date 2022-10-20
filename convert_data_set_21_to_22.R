library(data.table)
library(glue)

source("make_parameters.R")

files <- c("sim-donor_links3",
           "sim-donor_links4",
           "sim-Gsir_cov0_Dlidr-1",
           "sim-Gsir_cov0_Dlidr-2")

files <- c("fb-minmax")

for (f in files) {
    x <- readRDS(glue("param_sets/{f}.rds"))
    
    params <- make_parameters()
    
    for (n in names(params)) {
        if (n %in% names(x$params)) {
            params[[n]] <- x$params[[n]]
        }
    }
    
    params$sire_version <- "2.2"
    params$phi          <- 0.95
    
    if (is.null(x$protocol$seed)) {
        # No seed means testing coverage. Want 20 replicates and 4 threads each?
        # Compensate with fewer samples
        message("Setting up for coverage")
        params$nthreads <- 4L
        
        params$nsample      <- 2e6L
        params$burnin       <- 4e5L
        params$thin         <- 2e2L
    } else {
        # Fixed seed was for convergence. Only want 1st replicate and 10 threads
        message("Setting up for convergence")
        params$nthreads <- 10L
        x$protocol <- x$protocol[replicate == 1L]
        
        params$nsample      <- 1e7L
        params$burnin       <- 2e6L
        params$thin         <- 1e3L
    }
    
    if ("ignore_donors" %in% names(x$protocol)) {
        x$protocol[, ignore_donors := NULL]
    }
    
    x$params <- params
    
    x$protocol$data_set     <- paste0(x$protocol$data_set, "-mpi")
    x$protocol$sire_version <- x$params$sire_version
    x$protocol$nthreads     <- x$params$nthreads
    x$protocol$nsample      <- x$params$nsample
    x$protocol$burnin       <- x$params$burnin
    x$protocol$thin         <- x$params$thin
    x$protocol$phi          <- x$params$phi
    
    
    saveRDS(x, file = glue("param_sets/{f}-mpi.rds"))
}

for (f in files) {
    x <- readRDS(glue("param_sets/{f}-mpi.rds"))
    n <- nrow(x$protocol)
    t <- x$params$nthreads
    message(glue("{f} has {n} rows and {t} threads"))
}
