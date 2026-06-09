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

dataset <- "sim-test-inf2"

# Variable parameters ----
protocol <- rbind(
    # Misspecified models
    data.table(d = "FB_1_rpw, GEV SITTT, Weight SITTT, Fit d5, (Overfitting SITTT to none)"),    # 1
    data.table(d = "FB_1_rpw, GEV SITTT, Weight SITTT, Fit d6, (Overfitting SITTT to ST)"),      # 2
    data.table(d = "FB_1_rpw, GEV SITTT, Weight SITTT, Fit d7, (Underfitting SITTT to SILDT)"),  # 3
    data.table(d = "FB_1_rpw, GEV SITTT, Weight SITTT, Fit d8, (Underfitting SITTT to Cors=0)"), # 4
    data.table(d = "FB_1_rpw, GEV none,  Weight SITTT, Fit d1, (Underfitting none to SIT)"),     # 5
    data.table(d = "FB_1_rpw, GEV SITTT, Weight SITTT, Fit d1, (Testing h2=0)"),                 # 6

    fill = TRUE
)

# Setup
protocol[, setup := str_split_i(d, ", ", 1) |> str_to_lower(), .I]

# Genetic & Environmental Variance
protocol[, GEV := get_part(d, "GEV ") |> str_to_lower(), .I]
protocol[GEV == "sittt", `:=`(use_traits = "sildt", link_traits = "sittt")]
protocol[GEV != "sittt", use_traits := GEV]
protocol[, GEV := NULL]

protocol[, `:=`(vars = list(list(default = 0)),
                cors = list(list(default = 0)))]

# Handle FEs
protocol[, weight_fe := {
    x <- get_part(d, "Weight ")
    if (length(x) == 0 || x %in% c("", "none")) "" else str_to_lower(x)
}, .I]


# Handle patch_name: "Scen x" -> "scen-X-1"
protocol[, patch_name := get_part(d, "Fit ") |>
             str_replace("d(.*)", "scen-\\1-1"), .I]


# Common options ----
source("param_gen/common2.R")

common <- list(sim_new_data = "summary_sim",
               use_grm = "pedigree", # "HG_inv",
               inf_model = 4L,
               traits_source = "posterior", # should this be GRM?
               use_weight = "log",
               weight_is_nested = TRUE,
               cov_prior = list(type = "default", vals = c()),
               single_prior = "inverse",
               # expand_priors = 4,
               group_effect = 0.05,
               patch_dataset = "sim-test",
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

