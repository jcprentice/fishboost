library(data.table) # used for all data table manipulations
library(glue)

source("make_parameters.R")

# Overview ----

# Are we testing convergence or coverage?
convergence <- TRUE # only need convergence for FB data set
print(glue("convergence = {convergence}"))

data_set <- glue("sim-events-{x}",
                 x = if (convergence) 1 else 2)


# Variable parameters ----
protocol <- data.table(description = "FB Trial 1+2, SEIDR, parasites SEI, Traits SIR, FEs LIDR, GE 0.1",
                       setup = "fishboost",
                       use_parasites = "SEI",
                       censor = 0.8,
                       use_traits = "sir",
                       trial_fe = "lidr",
                       txd_fe = "lidr",
                       donor_fe = "lidr",
                       group_effect = 0.1,
                       prior__latent_period__type = "Fixed",
                       prior__latent_period__val1 = 10,
                       gamma_type = "gamma",
                       r_gamma_shape = 2,
                       time_step = 1,
                       use_fb_data = FALSE)


protocol[, label := paste0("s1", letters[1:.N])]

# Append "coverage" or "convergence" to description
if (TRUE) {
    co_str <- if (convergence) ", convergence" else ", coverage"
    protocol[, description := paste0(description, co_str)]
}

## Add replicates ----
n_replicates <- if (convergence) 1L else 20L
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
                          use_fb_data = FALSE)

params$use_fb_data <- FALSE
params$trial_fe <- "none"
params$donor_fe <- "none"
params$txd_fe <- "none"
params$group_effect <- -1
params$link_traits <- "slidr"
params$sim_link_trial <- "slidr"
params$sim_link_donor <- "slidr"
params$sim_link_shapes <- "ldr"
params$link_trial <- "slidr"
params$link_donor <- "slidr"
params$link_shapes <- "ldr"
params$pass_events <- 2
params$seed <- if (convergence) -1L else 1L
params$nthreads <- if (convergence) 16L else 4L
params$phi <- 0.95
nsample <- 1e6L
params$nsample <- as.integer(nsample)
params$burnin <- as.integer(nsample / 5)
params$thin <- as.integer(max(nsample / 1e4, 1))
params$nsample_per_gen <- as.integer(nsample * 2.5 / 1e3)
params$anneal <- "on"
params$anneal_power <- 4
params$sire_version <- "sire22"

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
                  "use_fb_data", "setup", "group_effect", "trial_fe", "donor_fe", "txd_fe",
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
