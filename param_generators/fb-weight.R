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

dataset <- "fb-weight"


# Variable parameters ----
protocol <- rbind(
    data.table(description = "FB_1, FEs ILDT"),
    data.table(description = "FB_2, FEs ILDT"),
    data.table(description = "FB_12, FEs ILDT"),
    data.table(description = "FB_1, FEs I[LDT]"),
    data.table(description = "FB_2, FEs I[LDT]"),
    data.table(description = "FB_12, FEs I[LDT]"),
    
    fill = TRUE
)

protocol[, setup := description |> str_split_i(", ", 1) |> str_to_lower()]

# Trial effects
protocol[setup == "fb_12", `:=`(trial_fe = "ildt",
                                txd_fe = "ildt")]

# Apply links
protocol[str_detect(description, "FEs ILDT"),
         `:=`(link_trial = "sildt", link_donor = "sildt", link_txd = "sildt", link_weight = "sildt")]
protocol[str_detect(description, "FEs I\\[LDT\\]"),
         `:=`(link_trial = "sittt", link_donor = "sittt", link_txd = "sittt", link_weight = "sittt")]

# Common options ----
source("param_generators/common2.R")

common <- list(use_traits = "sit",
               # trial_fe = "ildt",
               # txd_fe = "ildt",
               donor_fe = "ildt",
               weight_fe = "sildt",
               fix_donors = "no_Tsym_survivors",
               use_weight = "log",
               use_grm = "",
               group_effect = 0.1,
               prior__latent_period__type = "constant",
               prior__latent_period__val2 = 10,
               RP_dist = "gamma",
               time_step = 1) |>
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

