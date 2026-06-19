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

dataset <- "sim-events"

# Variable parameters ----
protocol <- rbind(
    data.table(d = "FB_1_rpw, Events 1, Fit d2", pass_events = "Tdeath"), # 1
    data.table(d = "FB_1_rpw, Events 2, Fit d2", pass_events = "Tsign,Tdeath"), # 2
    data.table(d = "FB_1_rpw, Events 3, Fit d2", pass_events = "Tinf,Tsign,Tdeath"), # 3
    data.table(d = "FB_1_rpw, Events 4, Fit d2", pass_events = "Tinf,Tinc,Tsign,Tdeath"), # 4

    fill = TRUE
)

# Setup
protocol[, setup := str_split_i(d, ", ", 1) |> str_to_lower(), .I]

# Common options ----
source("param_gen/common2.R")

common <- list(sim_new_data = "summary_sim",
               use_traits = "sildt",
               link_traits = "sittt",
               use_grm = "HG_inv", # "pedigree",
               inf_model = "S=pC",
               traits_source = "posterior", # should this be GRM?
               use_weight = "log",
               weight_fe = "sittt",
               weight_is_nested = TRUE,
               cov_prior = list(type = "default", vals = c()),
               single_prior = "inverse",
               # expand_priors = 4,
               group_effect = 0.05,
               patch_dataset = "sim-test",
               patch_name = "scen-2-1",
               patch_type = "sampled",
               bici_cmd = "inf",
               fix_donors = "no_Tsign_survivors",
               censor = 0.8,
               nsample = 2e6,
               nchains = 16,
               sample_states = 100,
               time_step_bici = 1) |>
    safe_merge(common2)

# Labels
protocol[, label := str_c("s", 1:.N)]

# Append "coverage" or "convergence" to description
protocol[, d := str_c(d, ", ", goal) |> str_squish()] |>
    setnames("d", "description")

## Add replicates ----
n_replicates <- 20
protocol[, scenario := .I]
protocol <- protocol[rep(1:.N, each = n_replicates)]
protocol[, replicate := 1:.N, scenario]
protocol[, dataset := dataset]
protocol[, name := str_c("scen-", scenario, "-", replicate)]
protocol[, seed := replicate]

# Set up patches
protocol[, patch_state := replicate]

# Prefer to have these columns in this order at the start
setcolorder(protocol, intersect(cols, names(protocol)))

message(str_glue("Protocol file '{dataset}' has:",
                 "- {nrow(protocol)} scenarios",
                 "- each with {ncol(protocol)} parameters",
                 "- and a further {length(common)} common parameters",
                 .sep = "\n"))

# Save to file ----
saveRDS(list(protocol = protocol,
             common = common),
        file = str_glue("param_sets/{dataset}.rds"))

