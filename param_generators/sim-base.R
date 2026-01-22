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

dataset <- "sim-base"

# Variable parameters ----
protocol <- rbind(
    data.table(d = "FB_1_rpw,  GEV SIT,   FE SIT, GRM HG_inv"), # 1
    data.table(d = "FB_12_rpw, GEV SITTT, FE SIT, GRM HG_inv"), # 2
    data.table(d = "FB_1_rpw,  GEV SITTT, FE SIT, GRM HG_inv"), # 2
    data.table(d = "FB_2_rpw,  GEV SILDT, FE SIT, GRM HG_inv"), # 6

    fill = TRUE
)

protocol[, `:=`(description = str_squish(d), d = NULL)]

protocol[, setup := description |> str_split_i(", ", 1), .I]

# Handle Genetic & Environmental Variance (GEV)
protocol[, GEV := description |> str_split_1(", ") |> str_subset("GEV") |>
             str_split_i(" ", 2) |> str_to_lower(), .I]
protocol[, use_traits := GEV]
protocol[GEV == "sittt", `:=`(use_traits = "sildt", link_traits = "sittt")]
protocol[, GEV := NULL]

# Handle GRM
protocol[, use_grm := description |> str_split_1(", ") |> str_subset("GRM") |>
             str_split_i(" ", 2), .I]
# protocol[, traits_source := fifelse(use_grm == "pedigree", "pedigree", "grm")]
protocol[, traits_source := "posterior"]

# Handle patch
protocol[, patch_name := description |> str_split_1(", ") |>
             str_subset("Fit to") |> str_replace_all("Fit to d(.*)", "scen-\\1-1"), .I]



protocol[, weight_fe := description |> str_split_1(", ") |> str_subset("FE") |>
             str_split_i(" ", 2) |> str_to_lower(), .I]

# Common options ----
source("param_generators/common2.R")

common <- list(sim_new_data = "bici",
               model_type = "SEIDR",
               use_grm = "pedigree",
               traits_source = "none",
               use_weight = "log",
               weight_is_nested = TRUE,
               # expand_priors = 4,
               group_effect = 0.05,
               patch_dataset = "fb-test",
               patch_name = "scen-1-1",
               patch_type = "median",
               patch_state = TRUE,
               # skip_patches = "beta", # "cov,base,beta",
               # cors <- c(si = -0.3, sl = 0.2, sd = 0.2, st = -0.3, il = 0.2,
               #           id =  0.2, it = 0.3, ld = 0.2, lt =  0.2, dt =  0.2),
               # latent_periods = 10,
               # detection_periods = 20,
               # removal_periods = 10,
               trial_fe = "sildt",
               donor_fe = "sildt",
               txd_fe = "sildt",
               bici_cmd = "sim",
               censor = 0.8,
               nsample = 1e4,
               sample_states = 100,
               nreps = 20,
               nchains = 1,
               time_step_bici = 0.2,
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

