{
    library(data.table)
    library(purrr)
    library(stringr)

    source("utils.R")
}

# Overview ----

# Are we testing convergence or coverage?
# (Coverage only makes sense with simulated data.)
n <- 1
goal <- c("convergence", "coverage")[[n]]
message("Goal = ", goal)

dataset <- str_glue("sim-events{x}",
                     x = if (goal == "convergence") 1 else 2)

# Variable parameters ----
protocol <- rbind(
    data.table(description = "Sim FB1, basic SEIDR, 4 events",
               setup = "fb_1",
               pass_events = "Tinf,Tlat,Tsym,Tdeath"),
    
    data.table(description = "Sim FB1, basic SEIDR, 3 events",
               setup = "fb_1",
               pass_events = "Tlat,Tsym,Tdeath"),
    
    data.table(description = "Sim FB1, basic SEIDR, 2 events",
               setup = "fb_1",
               pass_events = "Tsym,Tdeath"),
    
    data.table(description = "Sim FB1+2, basic SEIDR, 2 events",
               setup = "fb_12",
               pass_events = "Tsym,Tdeath"),
    
    fill = TRUE
)


# Labels
protocol[, label := str_c("s1", letters[1:.N])]


# Append "coverage" or "convergence" to description
protocol[, description := str_c(description, ", ", goal)]

## Add replicates ----
n_replicates <- if (goal == "convergence") 1L else 20L
protocol[, scenario := .I]
protocol <- protocol[rep(1:.N, each = n_replicates)]
protocol[, replicate := 1:.N, scenario]
protocol[, dataset := dataset]

# Fixed parameters ----

# Save params along with protocol so we know the defaults for all entries
params <- make_parameters(model_type = "SEIDR",
                          setup = "fb_1",
                          use_traits = "",
                          vars = 1.0,
                          cors = 0,
                          group_layout = "fishboost",
                          group_effect = -1,
                          sim_new_data = "r")

params$trial_fe <- ""
params$donor_fe <- ""
params$group_effect <- -1
params$link_traits <- "sittt"
params$sim_link_trial <- "sildt"
params$sim_link_donor <- "sildt"
params$sim_link_shapes <- "ldt"
params$link_trial <- "sildt"
params$link_donor <- "sildt"
params$link_shapes <- "ldt"
params$pass_events <- "Tsym,Tdeath"
params$seed <- if (goal == "convergence") 0 else -1
params$nchains <- if (goal == "convergence") 16 else 4
nsample <- 1e6L
params$nsample <- as.integer(nsample)
params$burnin <- as.integer(nsample / 5)
params$thin <- as.integer(max(nsample / 1e4, 1))
params$nsample_per_gen <- as.integer(max(nsample * 10 / 1e3))
params$anneal <- "on"
params$anneal_power <- 4
params$sire_version <- "sire22"

# Add missing columns ----

replace_NAs <- function(col) {
    # Make sure columns exist and overwrite NA values with something useful
    if (col %in% names(protocol)) {
        protocol[is.na(get(col)), (col) := params[[col]]]
    } else {
        protocol[, (col) := params[[col]]]
    }
}

missing_cols <- c("model_type", "dataset", "name", "use_traits", "vars", "cors",
                  "sim_new_data", "setup", "group_effect", "donor_fe", "trial_fe",
                  "weight_is_nested",
                  "seed", "sire_version", "nsample_per_gen", "anneal", "anneal_power",
                  "nsample", "burnin", "thin", "nchains")

walk(missing_cols, replace_NAs)

# Tidy up ----

# Give default values for data name, setup, cors, group_layout, ...
protocol[, name := str_c("scen-", scenario, "-", replicate)]

# Prefer to have these columns in this order at the start
cols <-  c("dataset", "description", "scenario", "replicate", "label",
           "name", "sim_new_data", "model_type", "setup")
setcolorder(protocol, intersect(cols, names(protocol)))

message(str_glue("protocol file '{dataset}' has {nrow(protocol)} rows x {ncol(protocol)} cols"))


# Save to file ----
saveRDS(list(protocol = protocol,
             params = params),
        file = str_glue("param_sets/{dataset}.rds"))

# fwrite(protocol, file = "protocol-sim.tsv", sep = "\t", quote = TRUE)
