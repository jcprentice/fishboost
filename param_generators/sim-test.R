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

dataset <- "sim-test"

# Variable parameters ----
protocol <- rbind(
    data.table(d = "FB_1_rpw,  GEV SIT,   FE SIT, Fit s1"),  # 1
    data.table(d = "FB_1_rpw,  GEV SITTT, FE SIT, Fit s2"),  # 2
    data.table(d = "FB_2_rpw,  GEV SITTT, FE SIT, Fit s6"),  # 3
    data.table(d = "FB_12_rpw, GEV SITTT, FE SIT, Fit s10"), # 4
    # Misspecified models
    data.table(d = "FB_1_rpw,  GEV none,  FE SIT, Fit s1"),  # 5
    data.table(d = "FB_1_rpw,  GEV ST,    FE SIT, Fit s1"),  # 6
    data.table(d = "FB_1_rpw,  GEV SILDT, FE SIT, Fit s1"),  # 7

    fill = TRUE
)

# Setup
protocol[, setup := str_split_i(d, ", ", 1) |> str_to_lower(), .I]

# Handle Genetic & Environmental Variance (GEV)
protocol[, GEV := get_part(d, "GEV") |> str_to_lower(), .I]
protocol[, use_traits := GEV]
protocol[GEV == "sittt", `:=`(use_traits = "sildt", link_traits = "sittt")]
protocol[, GEV := NULL]

# Handle fit
protocol[, patch_name := get_part(d, "Fit") |>
             str_replace("s(.*)", "scen-\\1-1"), .I]

# Handle GRM
protocol[, use_grm := get_part(d, "GRM"), .I]
# protocol[, traits_source := fifelse(use_grm == "pedigree", "pedigree", "grm")]
protocol[, traits_source := "posterior"]


protocol[, weight_fe := get_part(d, "FE") |> str_to_lower(), .I]

protocol[use_traits == "sildt", `:=`(vars = 0.5,
                                     cors = 0.2)]

# Skip patches when we want 0 variance
protocol[use_traits == "st", skip_patches := "cov_G_ii,cov_E_ii"]



# Common options ----
source("param_generators/common2.R")

common <- list(sim_new_data = "bici",
               model_type = "SEIDR",
               setup = "fb_12_rpw",
               use_grm = "posterior",
               inf_model = 4L,
               traits_source = "none",
               use_weight = "log",
               weight_is_nested = TRUE,
               # expand_priors = 4,
               group_effect = 0.05,
               patch_dataset = "fb-qtest",
               patch_type = "median",
               patch_state = FALSE,
               # skip_patches = "beta", # "cov,base,beta",
               # cors <- c(si = -0.3, sl = 0.2, sd = 0.2, st = -0.3, il = 0.2,
               #           id =  0.2, it = 0.3, ld = 0.2, lt =  0.2, dt =  0.2),
               # latent_periods = 10,
               # detection_periods = 20,
               # removal_periods = 10,
               trial_fe = "ildt",
               donor_fe = "ildt",
               txd_fe = "ildt",
               bici_cmd = "sim",
               censor = 0.8,
               nsample = 1e4,
               sample_states = 100,
               nreps = 20,
               time_step_bici = 0.2,
               ie_output = "true") |>
    safe_merge(common2)

common$nchains <- 1

# Labels
protocol[, label := str_c("s", 1:.N)]

# Append "coverage" or "convergence" to description
protocol[, d := str_c(str_squish(d), ", ", goal)] |>
    setnames("d", "description")


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

