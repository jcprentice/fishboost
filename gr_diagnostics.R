library(data.table)
library(gtools)
library(coda)
library(glue)
library(xtable)
library(stringr)

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
    
    cols <- parameter_estimates$parameter
    cols <- cols[!startsWith(cols, "Group effect")]
    
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

data_set <- "fb-mpi"
# data_set <- "sim-donor_links1-1-mpi"
scenario <- 4
replicate <- 0 # which replicate if doing coverage

xl <- list()
for (scenario in 1:4) {
    out <- gr_diagnostics(data_set, scenario, replicate)
    
    x <- as.data.table(round(out$GRD$psrf, 2), keep.rownames = "Parameter")
    x[, ESS := as.integer(out$ESS)]
    x[, Parameter := str_to_title(sub("_", " ", Parameter))]
    x[, Parameter := sub("Eta", "LP", Parameter)]
    x[, Parameter := sub("Rho", "DP", Parameter)]
    x[, Parameter := sub("Gamma", "RP", Parameter)]
    
    xl[[scenario]] <- x
    
    print(x)
    
    # print(xtable(x, align = "ll|rrr", caption = data_set), include.rownames = FALSE, booktabs = TRUE)
}
# 
