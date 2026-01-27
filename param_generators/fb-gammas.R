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

dataset <- "fb-gammas"


# Variable parameters ----
protocol <- rbind(
    # Standard
    data.table(d = "FB, Traits SIT,   FEs ILDT,               pedigree"), # 1
    data.table(d = "FB, Traits SITTT, FEs ILDT,               pedigree"), # 2
    # Replace var_E(tol)
    data.table(d = "FB, Traits SIT,   FEs ILDT, No VarE(T),   pedigree"), # 3
    data.table(d = "FB, Traits SITTT, FEs ILDT, No VarE(T),   pedigree"), # 4
    # Replace var_E(all)
    data.table(d = "FB, Traits SIT,   FEs ILDT, No VarE(SIT), pedigree"), # 5
    data.table(d = "FB, Traits SITTT, FEs ILDT, No VarE(SIT), pedigree"), # 6
    
    fill = TRUE
)

message(nrow(protocol), " scenarios")


# Set traits to SIT or SITTT
protocol[str_detect(d, "Traits SIT"),
         use_traits := "sit"]
protocol[str_detect(d, "Traits SITTT"),
         `:=`(use_traits = "all", link_traits = "sittt")]

protocol[str_detect(d, "No VarE\\(T\\)"),
         ge_opts := "no_ev_t"]

protocol[str_detect(d, "No VarE\\(SIT\\)"),
         ge_opts := "gt_only"]

# Common options ----
source("param_generators/common2.R")

common <- list(setup = "fb_12",
               trial_fe = "ildt",
               donor_fe = "ildt",
               txd_fe = "ildt",
               weight_fe = "sildt",
               weight_is_nested = TRUE,
               use_weight = "log",
               fix_donors = "no_Tsym_survivors",
               group_effect = 0.1,
               # expand_priors = 5,
               prior__latent_period__type = "constant",
               prior__latent_period__val2 = 5,
               prior__eta_shape__type = "constant",
               prior__eta_shape__val2 = 5,
               LP_type = "gamma",
               DP_dist = "gamma",
               RP_dist = "gamma",
               nsample = 1e6,
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

