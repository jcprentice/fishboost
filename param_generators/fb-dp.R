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

dataset <- "fb-dp"

# Variable parameters ----
protocol <- rbind(
    data.table(d = "FB_12_drop71, Traits SIT, FEs ILDT, DP 1,  GRM"), # 1
    data.table(d = "FB_12_drop71, Traits SIT, FEs ILDT, DP 10, GRM"), # 2
    data.table(d = "FB_12_drop71, Traits SIT, FEs ILDT, DP 20, GRM"), # 3
    data.table(d = "FB_12_drop71, Traits SIT, FEs ILDT, DP 30, GRM"), # 4
    data.table(d = "FB_12_drop71, Traits SIT, FEs ILDT, DP 40, GRM"), # 5
    data.table(d = "FB_12_drop71, Traits SIT, FEs ILDT, DP 50, GRM"), # 6
    
    fill = TRUE
)

message(nrow(protocol), " scenarios")

# Set traits to SIT or SITTT
protocol[str_detect(d, "Traits SIT"),
         use_traits := "sit"]
protocol[str_detect(d, "Traits SITTT"),
         `:=`(use_traits = "all", link_traits = "sittt")]

# Set LP
f <- function(x) x
protocol[, prior__detection_period__val2 := get_part(d, "DP") |> as.numeric(), .I]


# Set FEs to none, ILDT, or I[LDT]
protocol[str_detect(d, "FEs I\\[LDT\\]"),
         `:=`(link_trial = "sittt", link_donor = "sittt", link_txd = "sittt", link_weight = "sittt")]

# Common options ----
source("param_generators/common2.R")

common <- list(sim_new_data = "no",
               setup = "fb_12_drop71",
               use_grm = "H",
               prior__detection_period__type = "constant",
               fix_donors = "no_Tsym_survivors",
               t_demote = "10,80",
               trial_fe = "ildt",
               donor_fe = "ildt",
               txd_fe = "ildt",
               weight_fe = "sildt",
               weight_is_nested = TRUE,
               use_weight = "log",
               group_effect = 0.5,
               expand_priors = 4,
               nsample = 1e5,
               nsamples_per_gen = 1e5*3e-3,
               # sample_states = 100,
               ie_output = "true") |>
    safe_merge(common2)

# Labels
protocol[, label := str_c("s", 1:.N)]

# Append "coverage" or "convergence" to description
protocol[, d := str_c(d, ", ", goal) |> str_squish()] |>
    setnames("d", "description")

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

