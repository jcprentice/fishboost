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

dataset <- "sim-bici"

# Variable parameters ----
protocol <- rbind(
    data.table(description = "FB_12_rpw, Traits SIT, FEs ILDT, GRM pedigree, scen-5-1, timestep 0.05"), # 1
    data.table(description = "FB_12_rpw, Traits SIT, FEs ILDT, GRM pedigree, scen-5-1, timestep 0.1"), # 2
    data.table(description = "FB_12_rpw, Traits SIT, FEs ILDT, GRM pedigree, scen-5-1, timestep 0.2"), # 3
    data.table(description = "FB_12_rpw, Traits SIT, FEs ILDT, GRM pedigree, scen-5-1, timestep 0.5"), # 4
    
    fill = TRUE
)

protocol[, time_step_bici := description |> str_split_1(", ") |> str_subset("timestep") |>
             str_split_i(" ", 2) |> as.numeric(), .I]

# Common options ----
source("param_generators/common2.R")

common <- list(sim_new_data = "etc_sim",
               setup = "fb_12_rpw",
               use_traits = "sit",
               trial_fe = "ildt",
               donor_fe = "ildt",
               txd_fe = "ildt",
               weight_fe = "sildt",
               vars = 0.5,
               traits_source = "pedigree",
               use_weight = "log",
               weight_is_nested = TRUE,
               # expand_priors = 4,
               group_effect = 0.05,
               patch_dataset = "sim-base",
               patch_name = "scen-5-1",
               patch_type = "sampled",
               bici_cmd = "inf",
               censor = 0.8,
               nsample = 1e5,
               sample_states = 1e2,
               ie_output = "true") |>
    safe_merge(common2)

# Labels
protocol[, label := str_c("s", 1:.N)]

# Append "coverage" or "convergence" to description
protocol[, description := str_c(description, ", ", goal)]

## Add replicates ----
n_replicates <- if (goal == "convergence") 1 else 20
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

