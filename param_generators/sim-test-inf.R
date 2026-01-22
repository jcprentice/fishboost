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

dataset <- "sim-test-inf"

# Variable parameters ----
protocol <- rbind(
    data.table(d = "FB_1_rpw,  GEV SIT,   FEs SIT, Fit to d1, GRM HG_inv"), # 1
    data.table(d = "FB_12_rpw, GEV SITTT, FEs SIT, Fit to d2, GRM HG_inv"), # 2
    data.table(d = "FB_1_rpw,  GEV SITTT, FEs SIT, Fit to d3, GRM HG_inv"), # 3
    data.table(d = "FB_2_rpw,  GEV SITTT, FEs SIT, Fit to d4, GRM HG_inv"), # 4
    
    fill = TRUE
)

protocol[, `:=`(description = str_squish(d), d = NULL)]

protocol[, setup := str_split_i(description, ", ", 1) |> str_to_lower()]

# Handle patch_name: "Scen x" -> "scen-X-1"
protocol[, patch_name := description |> str_split_1(", ") |>
             str_subset("Fit to") |> str_replace_all("Fit to d(.*)", "scen-\\1-1"), .I]

# Common options ----
source("param_generators/common2.R")

common <- list(sim_new_data = "etc_sim",
               bici_cmd = "inf",
               use_grm = "HG_inv",
               use_traits = "sit",
               traits_source = "grm",
               cov_prior = "jeffreys",
               single_prior = "inverse",
               use_weight = "log",
               weight_is_nested = TRUE,
               # expand_priors = 4,
               patch_dataset = "sim-test",
               patch_type = "sampled",
               group_effect = 0.1,
               trial_fe = "ildt",
               donor_fe = "ildt",
               txd_fe = "ildt",
               weight_fe = "sit",
               censor = 0.8,
               nsample = 1e5,
               nchains = 8,
               sample_states = 100,
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

if (goal == "coverage") protocol[, seed := replicate]

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

