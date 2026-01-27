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

dataset <- "fb-ge"


# Variable parameters ----
protocol <- rbind(
    data.table(d = "FB_12, Traits SIT, FEs ILDT, GE 0.0, pedigree"),
    data.table(d = "FB_12, Traits SIT, FEs ILDT, GE 0.1, pedigree"),
    data.table(d = "FB_12, Traits SIT, FEs ILDT, GE 0.2, pedigree"),
    data.table(d = "FB_12, Traits SIT, FEs ILDT, GE 0.5, pedigree"),
    data.table(d = "FB_12, Traits SIT, FEs ILDT, GE 1.0, pedigree"),
    data.table(d = "FB_12, Traits SIT, FEs ILDT, GE 2.0, pedigree"),
    
    fill = TRUE
)

protocol[, group_effect := get_part(d, "GE") |> as.numeric(), .I]


# Common options ----
source("param_generators/common2.R")

common <- list(setup = "fb_12",
               use_traits = "sit",
               trial_fe = "ildt",
               donor_fe = "ildt",
               txd_fe = "ildt",
               weight_fe = "sildt",
               weight_is_nested = TRUE,
               use_weight = "log",
               fix_donors = "no_Tsym_survivors",
               prior__latent_period__type = "constant",
               prior__latent_period__val2 = 8,
               RP_dist = "gamma") |>
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

