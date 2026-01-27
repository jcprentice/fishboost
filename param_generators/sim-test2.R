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

dataset <- "sim-test2"

# Variable parameters ----
protocol <- rbind(
    data.table(d = "FB_12_rpw,    Traits SIT,  FEs ILDT, RP ~ exp, GRM HS_inv"), # 1
    data.table(d = "FB_12_rpw,    Traits SIT,  FEs ILDT, RP ~ exp, GRM HG_inv"), # 2
    data.table(d = "FB_12_drop71, Traits SIT,  FEs ILDT, RP ~ exp, GRM H_inv"), # 3
    data.table(d = "FB_12_drop71, Traits SIT,  FEs ILDT, RP ~ exp, GRM pedigree"), # 4
    data.table(d = "FB_12_drop71, Traits none,           RP ~ exp, GRM pedigree"), # 5

    fill = TRUE
)

protocol[, use_traits := get_part(d, "Traits") |> str_to_lower(), .I]

protocol[, setup := d |> str_split_i(", ", 1) |> str_to_lower(), .I]

protocol[, use_grm := get_part(d, "GRM"), .I]

protocol[str_detect(d, "pedigree"), traits_source := "pedigree"]
protocol[!str_detect(d, "pedigree"), traits_source := "grm"]

protocol[str_detect(d, "FEs ILDT"),
         `:=`(trial_fe = "ildt", donor_fe = "ildt", txd_fe = "ildt",
              weight_fe = "sildt")]
protocol[!str_detect(d, "FEs ILDT"),
         `:=`(trial_fe = "", donor_fe = "", txd_fe = "", weight_fe = "")]

# Common options ----
source("param_generators/common2.R")

common <- list(sim_new_data = "bici",
               vars = list(1, 1, 0.5),
               traits_source = "grm",
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
               sample_states = 1e3,
               ie_output = "true") |>
    safe_merge(common2)

# Labels
protocol[, label := str_c("s", 1:.N)]

# Append "coverage" or "convergence" to description
protocol[, d := str_c(d, ", ", goal) |> str_squish()] |>
    setnames("d", "description")

## Add replicates ----
n_replicates <- if (goal == "convergence") 5 else 20
protocol[, scenario := .I]
protocol <- protocol[rep(1:.N, each = n_replicates)]
protocol[, replicate := 1:.N, scenario]
protocol[, dataset := dataset]
protocol[, name := str_c("scen-", scenario, "-", replicate)]
protocol[, seed := replicate]

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

