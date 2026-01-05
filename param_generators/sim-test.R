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

dataset <- "sim-test"

# Variable parameters ----
protocol <- rbind(
    data.table(description = "FB_12_drop71, Traits SIT, FEs ILDT, RP ~ exp, pedigree"), # 1
    data.table(description = "FB_1_drop71, Traits SIT, FEs ILDT, RP ~ exp, pedigree"), # 6
    data.table(description = "FB_2_drop71, Traits SIT, FEs ILDT, RP ~ exp, pedigree"), # 7
    
    fill = TRUE
)

protocol[, setup := str_split_i(description, ", ", 1) |> str_to_lower()]

protocol[str_detect(description, "pedigree"), use_grm := ""]
protocol[str_detect(description, "GRM"), use_grm := "H"]

f <- function(s)  s |> str_split_1(", ") |> str_subset("Seed") |> str_split_i(" ", 2) |> as.numeric()
protocol[, seed := f(description), .I]

protocol[str_detect(description, "FEs ILDT"),
         `:=`(trial_fe = "ildt", donor_fe = "ildt", txd_fe = "ildt",
              weight_fe = "sildt")]

# Common options ----
source("param_generators/common2.R")

common <- list(sim_new_data = "r",
               vars = list(1,1,0.5),
               use_traits = "sit",
               use_weight = "log",
               weight_is_nested = TRUE,
               # expand_priors = 4,
               group_effect = 0.1,
               prior__latent_period__type = "Flat",
               prior__latent_period__val1 = 0,
               prior__latent_period__val2 = 10,
               prior__detection_period__val2 = 30,
               prior__trial_l__val1 = -8,
               prior__trial_l__val2 = 0,
               prior__trial_d__val1 = 0,
               prior__trial_d__val2 = 8,
               prior__donor_l__val1 = -2,
               prior__donor_l__val2 = 7,
               prior__txd_l__val1 = 0,
               prior__txd_l__val2 = 8,
               prior__txd_d__val1 = -4,
               prior__txd_d__val2 = 4,
               censor = 0.8,
               nsample = 1e5,
               # sample_states = 100,
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

