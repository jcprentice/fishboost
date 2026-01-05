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

dataset <- "fb-bici"

# Variable parameters ----
protocol <- rbind(
    data.table(description = "FB_12_drop71, Traits SIT, FEs ILDT, RP ~ exp, GRM H_inv"), # 1
    data.table(description = "FB_12_rpw, Traits SIT, FEs ILDT, RP ~ exp, GRM HS_inv"), # 2
    data.table(description = "FB_12_rpw, Traits SIT, FEs ILDT, RP ~ exp, GRM HG_inv"), # 3
    
    fill = TRUE
)

protocol[, setup := description |> str_split_1(", ") |> str_subset("FB") |>
             str_to_lower(), .I]

protocol[, use_grm := description |> str_split_1(", ") |> str_subset("GRM") |>
             str_split_i(" ", 2), .I]

protocol[str_detect(description, "Traits SIT"),
         use_traits := "sit"]
protocol[str_detect(description, "Traits SI\\[LDT\\]"),
         `:=`(use_traits = "sildt", link_traits = "sittt")]


# Common options ----
source("param_generators/common2.R")

common <- list(group_effect = 1.0,
               weight_is_nested = TRUE,
               use_weight = "log",
               # expand_priors = 4,
               trial_fe = "ildt",
               donor_fe = "ildt",
               txd_fe = "ildt",
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
               nsample = 1e5,
               sample_states = 1e2,
               ie_output = "true") |>
    safe_merge(common2)

# Labels
protocol[, label := str_c("s", 1:.N)]

# Append "coverage" or "convergence" to description
protocol[, description := str_c(description, ", ", goal)]

## Add replicates ----
n_replicates <- if (goal == "convergence") 10 else 20
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

