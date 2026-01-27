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

dataset <- "fb-final"

# Variable parameters ----
protocol <- rbind(
    # Pedigree
    data.table(d = "FB,       Traits SIT,   FEs ILDT, Weight nested, pedigree"), # 1
    data.table(d = "FB,       Traits SITTT, FEs ILDT, Weight nested, pedigree"), # 2
    data.table(d = "FB,       Traits SIT,   FEs ITTT, Weight nested, pedigree"), # 3
    data.table(d = "FB,       Traits SITTT, FEs ITTT, Weight nested, pedigree"), # 4
    # GRMs
    data.table(d = "FB,       Traits SIT,   FEs ILDT, Weight nested, GRM"), # 5
    data.table(d = "FB,       Traits SITTT, FEs ILDT, Weight nested, GRM"), # 6
    data.table(d = "FB,       Traits SIT,   FEs ITTT, Weight nested, GRM"), # 7
    data.table(d = "FB,       Traits SITTT, FEs ITTT, Weight nested, GRM"), # 8
    # Extra tests
    data.table(d = "FB, SIDR, Traits SIT,   FEs ILDT, Weight nested, GRM"), # 9
    data.table(d = "FB, SIDR, Traits SITTT, FEs ILDT, Weight nested, GRM"), # 10
    data.table(d = "FB, SIDR, Traits SIT,   FEs ITTT, Weight nested, GRM"), # 11
    data.table(d = "FB, SIDR, Traits SITTT, FEs ITTT, Weight nested, GRM"), # 12
    
    fill = TRUE
)

# Set model_type
protocol[, model_type := fifelse(str_detect(d, "SIDR"), "SIDR", "SEIDR")]

# Set traits to SIT or SITTT
protocol[str_detect(d, "Traits SIT"),
         use_traits := "sit"]
protocol[str_detect(d, "Traits SITTT"),
         `:=`(use_traits = "all", link_traits = "sittt")]

# Set FEs to ILDT or ITTT
protocol[str_detect(d, "FEs ITTT"),
         `:=`(link_trial = "sittt", link_donor = "sittt", link_txd = "sittt", link_weight = "sittt")]

# Set samples etc. for pedigree vs GRM
protocol[, nsample := 5e4]
protocol[str_detect(d, "GRM"),
         `:=`(use_grm = "Hinv", nsample = 0.4 * nsample)]
protocol[, nsample_per_gen := pmax(3e-3 * nsample, 1)]

# Fix Seeders
protocol[str_detect(d, "Fix Seeders"),
         `:=`(fix_donors = "time,no_Tsym_survivors", t_demote = 10, fix_eq_time = TRUE)]

# Common options ----
source("param_generators/common2.R")

common <- list(setup = "fb_12",
               trial_fe = "ildt",
               donor_fe = "ildt",
               txd_fe = "ildt",
               weight_fe = "sildt",
               weight_is_nested = TRUE,
               use_weight = "log",
               group_effect = 0.1,
               # expand_priors = 5,
               prior__latent_period__type = "constant",
               prior__latent_period__val2 = 5,
               RP_dist = "gamma",
               sample_states = 100,
               ie_output = "true",
               time_step = 1) |>
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

