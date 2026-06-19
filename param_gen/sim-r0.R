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

dataset <- "sim-r0"

# Variable parameters ----
protocol <- rbind(
    # Basic models
    data.table(d = "FB_1_rpw, GEV SIT, Weight SIT, Fit s1"),                 # 1
    data.table(d = "FB_1_rpw, GEV SIT, Weight SIT, Fit s1, Select sus 0.2"), # 2
    data.table(d = "FB_1_rpw, GEV SIT, Weight SIT, Fit s1, Select inf 0.2"), # 3
    data.table(d = "FB_1_rpw, GEV SIT, Weight SIT, Fit s1, Select tol 0.2"), # 4
    data.table(d = "FB_1_rpw, GEV SIT, Weight SIT, Fit s1, Select r0 0.2"),  # 5

    fill = TRUE
)

protocol[, `:=`(select_on = get_part(d, "Select", 2),
                select_top = get_part(d, "Select", 3) |> as.numeric()), .I]



# Common options ----
source("param_gen/common2.R")

common <- list(sim_new_data = "bici",
               setup = "fb_12_single",
               use_grm = "pedigree",
               trans_tree = "on",
               traits_source = "pedigree",
               model_type = "SEIDR",
               use_traits = "sit",
               link_traits = "sildt",
               inf_model = "S=pC",
               traits_source = "none",
               use_weight = "log",
               weight_fe = "sit",
               weight_is_nested = TRUE,
               group_effect = -1,
               patch_dataset = "fb-test",
               patch_name = "scen-1-1",
               patch_type = "median",
               patch_state = FALSE,
               I0 = 1,
               # prior__beta_Tr1__true_val = 3,
               # prior__beta_Tr2__true_val = 3,
               # skip_patches = "beta",
               trial_fe = "ildt",
               donor_fe = "ildt",
               txd_fe = "ildt",
               bici_cmd = "sim",
               censor = 0.8,
               nsample = 1e4,
               sample_states = 100,
               nreps = 50,
               time_step_bici = 0.2) |>
    safe_merge(common2)

common$nchains <- 1

# Labels
protocol[, label := str_c("s", 1:.N)]

# Append "coverage" or "convergence" to description
protocol[, d := str_c(str_squish(d), ", ", goal)] |>
    setnames("d", "description")


## Add replicates ----
n_replicates <- if (goal == "convergence") 1 else 50
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

