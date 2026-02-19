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
    # Basic models
    data.table(d = "FB_1_rpw,  GEV SIT,   FE SIT, Fit d1, GRM HG_inv"), # 1
    data.table(d = "FB_1_rpw,  GEV SITTT, FE SIT, Fit d2, GRM HG_inv"), # 2
    data.table(d = "FB_2_rpw,  GEV SITTT, FE SIT, Fit d3, GRM HG_inv"), # 3
    data.table(d = "FB_12_rpw, GEV SITTT, FE SIT, Fit d4, GRM HG_inv"), # 4
    # Misspecified models
    data.table(d = "FB_1_rpw,  GEV SIT,   FE SIT, Fit d2, GRM HG_inv, (Fitting SIT to SITTT)"),        # 5
    data.table(d = "FB_1_rpw,  GEV SITTT, FE SIT, Fit d1, GRM HG_inv, (Fitting SITTT to SIT)"),        # 6
    data.table(d = "FB_1_rpw,  GEV SIT,   FE SIT, Fit d5, GRM HG_inv, (Overfitting SIT to ST)"),       # 7
    data.table(d = "FB_1_rpw,  GEV SITTT, FE SIT, Fit d5, GRM HG_inv, (Overfitting SITTT to ST)"),     # 8
    data.table(d = "FB_1_rpw,  GEV SIT,   FE SIT, Fit d6, GRM HG_inv, (Underfitting SIT to SILDT)"),   # 9
    data.table(d = "FB_1_rpw,  GEV SITTT, FE SIT, Fit d6, GRM HG_inv, (Underfitting SITTT to SILDT)"), # 10

    fill = TRUE
)

# Setup
protocol[, setup := str_split_i(d, ", ", 1) |> str_to_lower(), .I]

# Genetic & Environmental Variance
protocol[, GEV := get_part(d, "GEV ") |> str_to_lower(), .I]
protocol[GEV == "sittt", `:=`(use_traits = "sildt", link_traits = "sittt")]
protocol[GEV != "sittt", use_traits := GEV]
protocol[, GEV := NULL]


# Handle FEs
protocol[, weight_fe := {
    x <- get_part(d, "FE ")
    if (length(x) == 0 || x %in% c("", "none")) "" else str_to_lower(x)
}, .I]


# Handle GRM
protocol[, use_grm := get_part(d, "GRM"), .I]
protocol[, traits_source := fifelse(use_grm == "pedigree", "pedigree", "grm")]


# Handle patch_name: "Scen x" -> "scen-X-1"
protocol[, patch_name := get_part(d, "Fit ") |>
             str_replace("d(.*)", "scen-\\1-1"), .I]


# Common options ----
source("param_generators/common2.R")

common <- list(sim_new_data = "etc_sim",
               traits_source = "posterior",
               use_weight = "log",
               weight_is_nested = TRUE,
               cov_prior = "jeffreys",
               vars = list(0),
               single_prior = "inverse",
               # expand_priors = 4,
               group_effect = 0.05,
               patch_dataset = "sim-test",
               patch_type = "sampled",
               bici_cmd = "inf",
               fix_donors = "no_Tsym_survivors",
               censor = 0.8,
               nsample = 1e6,
               nchains = 8,
               sample_states = 100,
               time_step_bici = 1,
               ie_output = "true") |>
    safe_merge(common2)

# Labels
protocol[, label := str_c("s", 1:.N)]

# Append "coverage" or "convergence" to description
protocol[, d := str_c(d, ", ", goal) |> str_squish()] |>
    setnames("d", "description")

## Add replicates ----
n_replicates <- 10
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

