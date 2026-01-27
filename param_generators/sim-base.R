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

dataset <- "sim-base"

# Variable parameters ----
protocol <- rbind(
    data.table(d = "FB_12_rpw, GEV SIT,   FE SIT,   beta med, GRM HG_inv"), # 1
    data.table(d = "FB_12_rpw, GEV SITTT, FE SIT,   beta med, GRM HG_inv"), # 2
    data.table(d = "FB_12_rpw, GEV SIT,   FE SIT,   beta q25, GRM HG_inv"), # 1
    data.table(d = "FB_12_rpw, GEV SITTT, FE SIT,   beta q25, GRM HG_inv"), # 2
    data.table(d = "FB_12_rpw, GEV ST,    FE SIT,   beta med, GRM HG_inv"), # 5
    data.table(d = "FB_12_rpw, GEV SILDT, FE SILDT, beta med, GRM HG_inv"), # 6
    
    fill = TRUE
)

# Handle Genetic & Environmental Variance (GEV)
protocol[, GEV := get_part(d, "GEV") |> str_to_lower(), .I]
protocol[, use_traits := GEV]
protocol[GEV == "sittt", `:=`(use_traits = "sildt", link_traits = "sittt")]
protocol[, GEV := NULL]

# Handle GRM
protocol[, use_grm := get_part(d, "GRM"), .I]
# protocol[, traits_source := fifelse(use_grm == "pedigree", "pedigree", "grm")]
protocol[, traits_source := "posterior"]


# Set beta to either q25 or median
protocol[str_detect(d, "beta q25"),
         `:=`(prior__beta_Tr1__true_val = 1.493645,
              prior__beta_Tr2__true_val = 1.719060)]
protocol[str_detect(d, "beta med"),
         `:=`(prior__beta_Tr1__true_val = 0.612829,
              prior__beta_Tr2__true_val = 0.889292)]

protocol[, weight_fe := get_part(d, "FE") |> str_to_lower(), .I]

protocol[, vars := fcase(
    use_traits == "st",    list(c(sus = 1.0, inf = 0.0, lat = 0.0, det = 0.0, tol = 0.5)),
    use_traits == "sit",   list(c(sus = 1.0, inf = 1.5, lat = 0.0, det = 0.0, tol = 0.5)),
    use_traits == "sildt", list(c(sus = 1.0, inf = 1.5, lat = 0.5, det = 0.5, tol = 0.5)),
    default = list(0))]

protocol[, skip_patches := fifelse(use_traits == "st", "beta,cov_G_ii,cov_E_ii", "beta")]


# Common options ----
source("param_generators/common2.R")

common <- list(sim_new_data = "bici",
               model_type = "SEIDR",
               setup = "fb_12_rpw",
               use_grm = "pedigree",
               traits_source = "none",
               use_weight = "log",
               weight_is_nested = TRUE,
               # expand_priors = 4,
               group_effect = 0.05,
               patch_dataset = "fb-test",
               patch_name = "scen-1-1",
               patch_type = "median",
               patch_state = TRUE,
               # skip_patches = "beta", # "cov,base,beta",
               # cors <- c(si = -0.3, sl = 0.2, sd = 0.2, st = -0.3, il = 0.2,
               #           id =  0.2, it = 0.3, ld = 0.2, lt =  0.2, dt =  0.2),
               # latent_periods = 10,
               # detection_periods = 20,
               # removal_periods = 10,
               trial_fe = "sildt",
               donor_fe = "sildt",
               txd_fe = "sildt",
               bici_cmd = "sim",
               censor = 0.8,
               nsample = 1e4,
               sample_states = 100,
               nreps = 20,
               time_step_bici = 0.2,
               ie_output = "true") |>
    safe_merge(common2)

common$nchains <- 1

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

