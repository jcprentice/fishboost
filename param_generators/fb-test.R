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

dataset <- "fb-test"

# Variable parameters ----
protocol <- rbind(
    data.table(d = "FB_12_rpw, GEV SIT,   GRM HG_inv"), # 1
    data.table(d = "FB_1_rpw,  GEV SIT,   GRM HG_inv"), # 2
    data.table(d = "FB_2_rpw,  GEV SIT,   GRM HG_inv"), # 3
    data.table(d = "FB_12_rpw, GEV SITTT, GRM HG_inv"), # 4
    data.table(d = "FB_1_rpw,  GEV SITTT, GRM HG_inv"), # 5
    data.table(d = "FB_2_rpw,  GEV SITTT, GRM HG_inv"), # 6
    
    fill = TRUE
)

protocol[, setup := get_part(d, "FB") |> str_to_lower(), .I]

protocol[str_detect(d, "GEV SIT"),
         use_traits := "sit"]
protocol[str_detect(d, "GEV SITTT"),
         `:=`(use_traits = "sildt", link_traits = "sittt")]

# protocol[, inf_model := get_part(d, "inf_model") |> as.integer(), .I]


# Common options ----
source("param_generators/common2.R")

common <- list(use_grm = "HG_inv",
               popn_format = "intervals",
               inf_model = 1L,
               group_effect = 0.05,
               weight_is_nested = TRUE,
               cov_prior = "jeffreys",
               use_weight = "log",
               # expand_priors = 4,
               trial_fe = "ildt",
               donor_fe = "ildt",
               txd_fe = "ildt",
               weight_fe = "sit",
               prior__beta_Tr1__val2 = 4,
               prior__beta_Tr2__val2 = 3,
               prior__weight2_i__val1 = -2,
               prior__weight2_i__val2 = +6,
               prior__weight2_l__val1 = -2,
               prior__weight2_l__val2 = +6,
               `prior__latent_period_Tr2,Don__val1` = 1,
               `prior__latent_period_Tr2,Don__val2` = 20,
               fix_donors = "no_Tsym_survivors",
               nsample = if (str_detect(dataset, "q")) 5e5 else 1e7,
               sample_states = if (str_detect(dataset, "q")) 1e2 else 1e3,
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

