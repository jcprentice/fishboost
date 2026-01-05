{
    library(data.table)
    library(purrr)
    library(stringr)

    source("utils.R")
}

# Overview ----

# Are we testing convergence or coverage?
# (Coverage only makes sense with simulated data.)
n <- 2
goal <- c("convergence", "coverage")[[n]]
message("Goal = ", goal)

dataset <- str_c("sim-simple", n)


# Variable parameters ----
protocol <- rbind(
    # Pass all events
    data.table(description = "FB Trial 1+2, No IEs, No FEs, Pass all events"), # 1
    data.table(description = "FB Trial 1+2, FEs ILDT, Weight Nested SILDT, Pass all events"), # 2
    data.table(description = "FB Trial 1+2, IEs ILDT, No FEs, Pass all events"), # 3
    data.table(description = "FB Trial 1+2, IEs ILDT, FEs ILDT, Weight Nested SILDT, Pass all events"), # 4
    # Pass regular events
    data.table(description = "FB Trial 1+2, No IEs, No FEs"), # 5
    data.table(description = "FB Trial 1+2, FEs ILDT, Weight Nested SILDT"), # 6
    data.table(description = "FB Trial 1+2, IEs ILDT, No FEs"), # 7
    data.table(description = "FB Trial 1+2, IEs ILDT, FEs ILDT, Weight Nested SILDT"), # 12
    
    fill = TRUE
)

protocol[, use_traits := fifelse(str_detect(description, "IEs SIT"), "sit", "")]
protocol[str_detect(description, "FEs SIT"),
         `:=`(trial_fe = "ildt", donor_fe = "ildt", txd_fe = "ildt", weight_fe = "sildt")]
protocol[, pass_events := fifelse(str_detect(description, "Pass all events"), "all", NA)]


# Common options
protocol[, `:=`(sim_new_data = "r",
                patch_dataset = "",
                patch_type = "mean",
                traits_source = "pedigree",
                setup = "fb_12",
                weight_is_nested = TRUE,
                use_weight = "log",
                # use_grm = "",
                group_effect = 0.1,
                prior__latent_period__type = "Fixed",
                prior__latent_period__val2 = 10,
                expand_priors = 5,
                RP_dist = "gamma",
                sample_states = fifelse(goal == "convergence", 100, 0),
                ie_output = "false",
                time_step = 1)]

# expand_priors must be moved to before any other priors, or be overwritten
if ("expand_priors" %in% names(protocol)) {
    setcolorder(protocol, "expand_priors", before = str_which(names(protocol), "prior__")[[1]])
}

# Labels
protocol[, label := str_c("s", 1:.N)]


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
                          setup = "fb_12",
                          use_traits = "sit",
                          vars = 1.0,
                          cors = 0,
                          group_layout = "fishboost",
                          group_effect = -1,
                          sim_new_data = "r")

params$trial_fe <- ""
params$donor_fe <- ""
params$txd_fe <- ""
params$group_effect <- -1
params$link_traits <- "sildt"
params$sim_link_trial <- "sildt"
params$sim_link_donor <- "sildt"
params$sim_link_shapes <- "ldt"
params$link_trial <- "sildt"
params$link_donor <- "sildt"
params$link_shapes <- "ldt"
params$pass_events <- "Tsym,Tdeath"
params$seed <- if (goal == "convergence") 0 else -1
params$nchains <- if (goal == "convergence") 16 else 4
params$phi <- 1.0
nsample <- 1e6L
params$nsample <- as.integer(nsample)
params$burnin <- as.integer(nsample / 5)
params$thin <- as.integer(max(nsample / 1e4, 1))
params$nsample_per_gen <- as.integer(max(nsample * 2.5e-3, 1))
params$anneal <- "on"
params$anneal_power <- 4

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
                  "sim_new_data", "setup", "group_effect", "trial_fe", "donor_fe", "txd_fe",
                  "weight_is_nested",
                  "seed", "nsample_per_gen", #"anneal", "anneal_power",
                  "phi", "nsample", "burnin", "thin", "nchains")

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
