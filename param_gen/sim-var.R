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

dataset <- "sim-var"

# Variable parameters ----
protocol <- rbind(
    # Basic models
    data.table(d = "FB_1_rpw, GEV none, Cors 0"), # 1
    data.table(d = "FB_1_rpw, GEV S,    Cors 0"), # 2
    data.table(d = "FB_1_rpw, GEV I,    Cors 0"), # 3
    data.table(d = "FB_1_rpw, GEV T,    Cors 0"), # 4
    data.table(d = "FB_1_rpw, GEV SI,   Cors 0"), # 5
    data.table(d = "FB_1_rpw, GEV ST,   Cors 0"), # 6
    data.table(d = "FB_1_rpw, GEV IT,   Cors 0"), # 7
    data.table(d = "FB_1_rpw, GEV SIT,  Cors 0"), # 8

    data.table(d = "FB_1_rpw, GEV none, Cors 0.3"), # 9
    data.table(d = "FB_1_rpw, GEV S,    Cors 0.3"), # 10
    data.table(d = "FB_1_rpw, GEV I,    Cors 0.3"), # 11
    data.table(d = "FB_1_rpw, GEV T,    Cors 0.3"), # 12
    data.table(d = "FB_1_rpw, GEV SI,   Cors 0.3"), # 13
    data.table(d = "FB_1_rpw, GEV ST,   Cors 0.3"), # 14
    data.table(d = "FB_1_rpw, GEV IT,   Cors 0.3"), # 15
    data.table(d = "FB_1_rpw, GEV SIT,  Cors 0.3"), # 16

    fill = TRUE
)

# Setup
protocol[, setup := str_split_i(d, ", ", 1) |> str_to_lower(), .I]

# Handle Genetic & Environmental Variance (GEV)
protocol[, GEV := get_part(d, "GEV") |> str_to_lower(), .I]
protocol[, use_traits := GEV]
protocol[GEV == "sittt", `:=`(use_traits = "sildt", link_traits = "sittt")]
protocol[, GEV := NULL]

# Handle GRM
protocol[, traits_source := "pedigree"]


# Vars and Cors
protocol[, `:=`(vars = list(list(default = 0.5)),
                cors = list(list(default = 0)))]
protocol[str_detect(d, "Cors 0.3"), cors := list(list(default = 0.3))]



# Common options ----
source("param_gen/common2.R")

common <- list(sim_new_data = "bici",
               model_type = "SEIDR",
               use_grm = "pedigree", # "HG_inv",
               # setup = "fb_12_rpw",
               inf_model = 4L,
               traits_source = "none",
               use_weight = "log",
               weight_fe = "sit",
               weight_is_nested = TRUE,
               # expand_priors = 4,
               group_effect = 0.05,
               patch_dataset = "fb-test",
               patch_name = "scen-1-1",
               patch_type = "median",
               patch_state = FALSE,
               skip_patches = "cov",
               # skip_patches = "beta", # "cov,base,beta",
               # cors <- c(si = -0.3, sl = 0.2, sd = 0.2, st = -0.3, il = 0.2,
               #           id =  0.2, it = 0.3, ld = 0.2, lt =  0.2, dt =  0.2),
               # LP = 10,
               # DP = 20,
               # RP = 10,
               trial_fe = "ildt",
               donor_fe = "ildt",
               txd_fe = "ildt",
               bici_cmd = "sim",
               censor = 0.8,
               nsample = 1e4,
               sample_states = 100,
               nreps = 20,
               time_step_bici = 0.2) |>
    safe_merge(common2)

common$nchains <- 1

# Labels
protocol[, label := str_c("s", 1:.N)]

# Append "coverage" or "convergence" to description
protocol[, d := str_c(str_squish(d), ", ", goal)] |>
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
