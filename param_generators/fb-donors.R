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

dataset <- "fb-donors"

# Variable parameters ----
protocol <- rbind(
    # No donor fixes
    data.table(d = "FB_12_drop71, Traits SIT, FEs ILDT, RP ~ exp, leave donors, GRM"), # 1
    data.table(d = "FB_1_drop71,  Traits SIT, FEs ILDT, RP ~ exp, leave donors, GRM"), # 2
    data.table(d = "FB_2_drop71,  Traits SIT, FEs ILDT, RP ~ exp, leave donors, GRM"), # 3
    # Standard Fix
    data.table(d = "FB_12_drop71, Traits SIT, FEs ILDT, RP ~ exp, leave donors, GRM"), # 4
    data.table(d = "FB_1_drop71,  Traits SIT, FEs ILDT, RP ~ exp, leave donors, GRM"), # 5
    data.table(d = "FB_2_drop71,  Traits SIT, FEs ILDT, RP ~ exp, leave donors, GRM"), # 6
    # Strict donor fixes
    data.table(d = "FB_12_drop71, Traits SIT, FEs ILDT, RP ~ exp, strict donors, GRM"), # 7
    data.table(d = "FB_1_drop71,  Traits SIT, FEs ILDT, RP ~ exp, strict donors, GRM"), # 8
    data.table(d = "FB_2_drop71,  Traits SIT, FEs ILDT, RP ~ exp, strict donors, GRM"), # 9
    
    fill = TRUE
)

protocol[, setup := str_split_i(d, ", ", 1) |> str_to_lower()]

protocol[str_detect(d, "pedigree"), use_grm := ""]
protocol[str_detect(d, "GRM"), use_grm := "H"]

protocol[str_detect(d, "FEs ILDT"),
         `:=`(trial_fe = "ildt", donor_fe = "ildt", txd_fe = "ildt",
              weight_fe = "sildt")]

# All should have no Tsym survivors, but also check more stringent version
protocol[, fix_donors := fcase(str_detect(d, "leave donors"), "",
                               str_detect(d, "strict donors"), "no_Tsym_survivors,time",
                               default = "no_Tsym_survivors")]

# Common options ----
source("param_generators/common2.R")

common <- list(use_traits = "sit",
               group_effect = 1.0,
               weight_is_nested = TRUE,
               use_weight = "log",
               # expand_priors = 4,
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
               nsample = 1e4,
               # sample_states = 100,
               ie_output = "true") |>
    safe_merge(common2)

# Labels
protocol[, label := str_c("s", 1:.N)]

# Append "coverage" or "convergence" to description
protocol[, d := str_c(d, ", ", goal) |> str_squish()] |>
    setnames("d", "description")

## Add replicates ----
n_replicates <- if (goal == "convergence") 3 else 20
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

