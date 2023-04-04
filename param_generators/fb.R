library(data.table) # used for all data table manipulations
library(glue)

source("make_parameters.R")

# Overview ----

# Are we testing convergence or coverage?
convergence <- T
# mpi = SIRE 2.2, otherwise SIRE 2.1
sire22 <- T
print(glue("convergence = {convergence}, sire22 = {sire22}"))

data_set <- glue("fb{x}", x = if (sire22) "-mpi" else "")

# Variable parameters ----
protocol <- rbind(
    data.table(description = "FB Trial 1, SEIDR, Traits SIR, Donor LIDR, GE 0.3",
               setup = "fb1",
               label = "s1a",
               use_traits = "sir",
               link_traits = "slidr",
               link_donor = "slidr",
               donor_fe = "lidr",
               group_effect = 0.2),
    
    data.table(description = "FB Trial 2, SEIDR, Traits SIR, Donor LIDR, GE 0.3",
               setup = "fb2",
               label = "s1a",
               use_traits = "sir",
               link_traits = "slidr",
               link_donor = "slidr",
               donor_fe = "lidr",
               group_effect = 0.2),
    
    data.table(description = "FB Trial 1+2, SEIDR, Traits SIR, Donor LIDR, Trial LIDR, GE 0.3",
               setup = "fishboost",
               label = "s1c",
               use_traits = "sir",
               link_traits = "slidr",
               link_trial = "slidr",
               link_donor = "slidr",
               trial_fe = "lidr",
               donor_fe = "lidr",
               group_effect = 0.2),
    
    data.table(description = "FB Trial 1, basic",
               setup = "fishboost",
               label = "s1d",
               use_traits = "none",
               link_traits = "slidr",
               link_trial = "slidr",
               link_donor = "slidr",
               trial_fe = "none",
               donor_fe = "none",
               group_effect = -1),
    
    
    
    fill = TRUE
)

# Append "coverage" or "convergence" to description
if (TRUE) {
    co_str <- if (convergence) ", convergence" else ", coverage"
    protocol[, description := paste0(description, co_str)]
}

## Add replicates ----
n_replicates <- if (convergence) {if (sire22) 1L else 10L} else 20L
protocol[, scenario := .I]
protocol <- protocol[rep(1:.N, each = n_replicates)]
protocol[, replicate := 1:n_replicates, by = scenario]
protocol[, data_set := data_set]

# Fixed parameters ----

# Save params along with protocol so we know the defaults for all entries
params <- make_parameters(model_type = "SEIDR",
                          setup = "fishboost",
                          use_traits = "sir",
                          vars = 1.0,
                          covars = "none",
                          group_layout = "fishboost",
                          group_effect = -1,
                          use_fb_data = TRUE)

params$use_fb_data <- TRUE
params$trial_fe <- "none"
params$donor_fe <- "none"
params$group_effect <- -1
params$link_traits <- "srirr"
params$sim_link_trial <- "slidr"
params$sim_link_donor <- "slidr"
params$sim_link_shapes <- "ldr"
params$link_trial <- "slidr"
params$link_donor <- "slidr"
params$link_shapes <- "ldr"
params$pass_events <- 2
params$seed <- if (convergence) -1L else 1L
params$nthreads <- if (convergence) 10L else 4L
params$phi <- 0.95
nsample <- 1e6L
params$nsample <- as.integer(nsample)
params$burnin <- as.integer(nsample / 5)
params$thin <- as.integer(max(nsample / 1e5, 1)) # only fb wants 1e5 samples
params$nsample_per_gen <- as.integer(nsample * 3 / 1e3)
params$anneal <- "on"
params$anneal_power <- 4
params$sire_version <- if (sire22) "sire22" else "sire21"

# Add missing columns ----

replace_NAs <- function(col, val = "") {
    # Make sure columns exist and overwrite NA values with something useful
    if (col %in% names(protocol)) {
        protocol[is.na(get(col)), (col) := val]
    } else {
        protocol[, (col) := val]
    }
}

missing_cols <- c("model_type", "data_set", "name", "use_traits", "vars", "covars",
                  "use_fb_data", "setup", "group_effect", "donor_fe", "trial_fe", "txd_fe",
                  "seed", "sire_version", "nsample_per_gen", #"anneal", "anneal_power",
                  "phi", "nsample", "burnin", "thin", "nthreads")
for (mc in missing_cols) {
    replace_NAs(mc, params[[mc]])
}

# Tidy up ----

# Give default values for data name, setup, covars, group_layout, ...
protocol[, name := paste0("scen-", scenario, "-", replicate)]

# Prefer to have these columns in this order at the start
cols <-  c("data_set", "description", "scenario", "replicate", "label",
           "name", "use_fb_data", "model_type", "setup")
setcolorder(protocol, intersect(cols, names(protocol)))

message(glue("protocol file '{data_set}' has {nrow(protocol)} rows x {ncol(protocol)} cols"))


# Save to file ----
saveRDS(list(protocol = protocol,
             params = params),
        file = glue("param_sets/{data_set}.rds"))

# fwrite(protocol, file = "protocol-sim.tsv", sep = "\t", quote = TRUE)
