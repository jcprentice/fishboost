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

dataset <- "fb-fes3"


# Variable parameters ----
protocol <- rbind(
    data.table(d = "FB_1,  FEs ILDT, Weight"), # 1
    data.table(d = "FB_2,  FEs ILDT, Weight"), # 2
    data.table(d = "FB_12, FEs ILDT, Weight"), # 3
    data.table(d = "FB_12, FEs ILDT, Weight nested"), # 4
    data.table(d = "FB_1,  FEs IDT,  Weight"), # 5
    data.table(d = "FB_2,  FEs IDT,  Weight"), # 6
    data.table(d = "FB_12, FEs IDT,  Weight"), # 7
    data.table(d = "FB_12, FEs IDT,  Weight nested"), # 8
    data.table(d = "FB_1,  FEs IDT,  No weight"), # 9
    data.table(d = "FB_2,  FEs IDT,  No weight"), # 10
    data.table(d = "FB_12, FEs IDT,  No weight"), # 11
    data.table(d = "FB_12, FEs ILDT, Weight nested, GRM"), # 12
    data.table(d = "FB_12, FEs ILDT, Weight, GRM"), # 13
    
    fill = TRUE
)

protocol[, setup := d |> str_split_i(", ", 1) |> str_to_lower()]
protocol[, donor_fe := str_part(d, "FEs") |> str_to_lower()]
protocol[, `:=`(trial_fe = fifelse(str_detect(setup, "12"), donor_fe, ""),
                txd_fe = trial_fe,
                weight_fe = str_c("s", donor_fe))]
protocol[, weight_is_nested := str_detect(d, "Weight nested")]

# Common options ----
source("param_generators/common2.R")

common <- list(use_traits = "sit",
               fix_donors = "no_Tsym_survivors",
               use_weight = "log",
               # use_grm = "",
               group_effect = 0.1,
               prior__latent_period__type = "constant",
               prior__latent_period__val2 = 10,
               expand_priors = 5,
               RP_dist = "gamma",
               nsample = 1e6,
               burnprop = 0.2,
               sample_states = 100,
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

