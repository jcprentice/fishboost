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

dataset <- "fb-qtest"

# Variable parameters ----
protocol <- rbind(
    data.table(d = "FB_12_rpw, GEV SIT,   FE SIT, LP_D2 low,  GRM HG_inv"), # 1
    data.table(d = "FB_12_rpw, GEV SIT,   FE SIT, LP_D2 high, GRM HG_inv"), # 2
    data.table(d = "FB_12_rpw, GEV SITTT, FE SIT, LP_D2 low,  GRM HG_inv"), # 3
    data.table(d = "FB_12_rpw, GEV SITTT, FE SIT, LP_D2 high, GRM HG_inv"), # 4
    
    fill = TRUE
)

protocol[, `:=`(description = str_squish(d), d = NULL)]

protocol[str_detect(description, "GEV SIT"),
         use_traits := "sit"]
protocol[str_detect(description, "GEV SITTT"),
         `:=`(use_traits = "sildt", link_traits = "sittt")]

protocol[, weight_fe := description |> str_split_1(", ") |> str_subset("FE") |>
             str_split_i(" ", 2) |> str_to_lower(), .I]

protocol[str_detect(description, "LP_D2 low"),
         `:=`(`prior__latent_period_Tr2,Don__val1` = 1,
              `prior__latent_period_Tr2,Don__val2` = 20)]

protocol[str_detect(description, "LP_D2 high"),
         `:=`(`prior__latent_period_Tr2,Don__val1` = 20,
              `prior__latent_period_Tr2,Don__val2` = 120)]

# Common options ----
source("param_generators/common2.R")

common <- list(setup = "fb_12_rpw",
               use_grm = "HG_inv",
               popn_format = "intervals",
               group_effect = 0.05,
               weight_is_nested = TRUE,
               cov_prior = "jeffreys",
               use_weight = "log",
               # expand_priors = 4,
               trial_fe = "ildt",
               donor_fe = "ildt",
               txd_fe = "ildt",
               prior__beta_Tr1__val2 = 1,
               prior__beta_Tr2__val2 = 1,
               prior__beta_Tr1__val2 = 1,
               prior__cov_G_ii__val2 = 4,
               prior__weight2_i__val1 = -2,
               prior__weight2_i__val2 = 6,
               prior__weight2_l__val1 = -2,
               prior__weight2_l__val2 = 6,
               fix_donors = "no_Tsym_survivors",
               nsample = 2e5,
               sample_states = 1e2,
               ie_output = "true") |>
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

