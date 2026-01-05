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

dataset <- "fb-weight-grm"


# Variable parameters ----
protocol <- rbind(
    data.table(description = "FB_12, GRM",
               use_grm = "Hinv"),
    
    data.table(description = "FB_12, pedigree",
               use_grm = ""),
    
    fill = TRUE
)


# Common options ----
source("param_generators/common2.R")

common <- list(setup = "fb_12",
               use_traits = "sit",
               trial_fe = "ildt",
               txd_fe = "ildt",
               donor_fe = "ildt",
               weight_fe = "sildt",
               link_trial = "sildt",
               link_donor = "sildt",
               link_txd = "sildt",
               link_weight = "sildt",
               fix_donors = "no_Tsym_survivors",
               use_weight = "log",
               group_effect = 0.1,
               prior__latent_period__type = "constant",
               prior__latent_period__val2 = 10,
               RP_dist = "gamma") |>
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

