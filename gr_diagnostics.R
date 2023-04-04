{
    library(data.table)
    library(gtools)
    library(coda)
    library(glue)
    library(xtable)
    library(stringr)
    
    source("rename_pars.R")
}

gr_diagnostics <- function(data_set = "sim", scenario = 1, replicate = 0) {
    
    # data_set = "sim-events1"; scenario = 1; replicate = 1
    
    data_dir <- glue("data/{data_set}")
    res_dir  <- glue("results/{data_set}")

    trace_files <- list.files(path = data_dir,
                              pattern = "trace",
                              recursive = TRUE,
                              full.names = TRUE) |>
        mixedsort(decreasing = TRUE)
    
    # Sieve for the right scenario
    trace_files <- trace_files[grepl(glue("scen-{scenario}-"), trace_files)]
    
    # Drop combine files for SIRE 2.2
    trace_files <- trace_files[!grepl("combine", trace_files)]
    
    # If doing coverage, then maybe specify which replicate to keep
    if (replicate > 0) {
        trace_files <- trace_files[grepl(glue("-{replicate}_out"), trace_files)]
    }
    
    res <- list.files(res_dir, pattern = glue("scen-{scenario}-"), full.names = TRUE)
    load(res[1])
    
    # Select columns to keep
    cols <- parameter_estimates$parameter
    cols <- cols[!startsWith(cols, "Group effect")]
    # cols <- cols[!startsWith(cols, "L_")]
    # cols <- cols[!cols %in% c("state", "Prior", "Posterior", "Number infected", "log(phi)")]
    
    # Remove columns with zero variance
    cols <- cols[!cols %in% params$priors[type == "Fixed", parameter]]
    
    
    thin <- params$thin
    burnin <- params$burnin
    bstart <- burnin / thin + 1
    
    y <- data.table()
    l <- mcmc.list()
    for (i in seq_along(trace_files)) {
        x <- fread(trace_files[i])
        
        l[[i]] <- mcmc(x[, ..cols], start = bstart, thin = thin)
        y <- rbind(y, x[bstart:.N, ..cols])
    }
    
    list(ESS = effectiveSize(y),
         GRD = gelman.diag(l))
}

# data_set <- "sim-events1-mpi"; nScenarios <- 4
# data_set <- "sim-donor_links1-1"; nScenarios <- 3
# data_set <- "sim-donor_links1-1-mpi"; nScenarios <- 4
# data_set <- "sim-G_Da-1-mpi"; nScenarios <- 4
# data_set <- "sim-Gsi_cov_Da-1-mpi"; nScenarios <- 4
# data_set <- "fb-mpi"; nScenarios <- 4
data_set <- "fb-parasites3"; nScenarios <- 12

replicate <- 0 # which replicate if doing coverage

xl <- list()
gd <- list()
for (scenario in 1:nScenarios) {
    out <- gr_diagnostics(data_set, scenario, replicate)
    gd[[scenario]] <- out$GRD$mpsrf
    
    x <- as.data.table(round(out$GRD$psrf, 2), keep.rownames = "Parameter")
    x[, `Point est.` := NULL]
    x[, ESS := format(as.integer(out$ESS), big.mark = ",")]
    x[, Parameter := rename_pars(Parameter)]
    setnames(x, "Upper C.I.", "GR95")
    
    x <- rbind(x, data.table(Parameter = "MVPSF", GR95 = gd[[scenario]]), fill = TRUE)
    xl[[scenario]] <- x
    
    print(x)
}

print(round(unlist(gd), 2))

trial_offset <- 0:2
names(trial_offset) <- c("1", "2", "12")

{
    offset <- trial_offset[["12"]]
    foo <- cbind(xl[[1 + offset]][, .(Parameter)],
                 xl[[1 + offset]][, .(GR95, ESS)],
                 xl[[4 + offset]][, .(GR95, ESS)],
                 xl[[7 + offset]][, .(GR95, ESS)],
                 xl[[10 + offset]][, .(GR95, ESS)])
    
    names(foo) <- c("Parameter", rep(c("95% GR", "ESS"), 4))
    
    print(xtable(foo, align = "llrrrrrrrr"),
          include.rownames = FALSE,
          booktabs = TRUE)
}
