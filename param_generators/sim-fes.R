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

dataset <- "sim-fes"


# Variable parameters ----
protocol <- rbind(
    data.table(description = "FB Trial 1, FEs IDR, Weight", # 5
               patch_name = "scen-5-1",
               setup = "fb_1",
               donor_fe = "idr",
               weight_fe = "sidr"),
    
    data.table(description = "FB Trial 2, FEs IDR, Weight", # 6
               patch_name = "scen-6-1",
               setup = "fb_2",
               donor_fe = "idr",
               weight_fe = "sidr"),
    
    data.table(description = "FB Trial 1+2, FEs IDR, Weight", # 7
               patch_name = "scen-7-1",
               setup = "fb_12",
               donor_fe = "idr",
               trial_fe = "idr",
               txd_fe = "idr",
               weight_fe = "sidr",
               weight_is_nested = FALSE),
    
    data.table(description = "FB Trial 1+2, FEs IDR, Weight nested", # 8
               patch_name = "scen-8-1",
               setup = "fb_12",
               donor_fe = "idr",
               trial_fe = "idr",
               txd_fe = "idr",
               weight_fe = "sidr",
               weight_is_nested = TRUE),
    
    data.table(description = "FB Trial 1+2, FEs IDR, Weight nested", # 8
               patch_name = "scen-8-1",
               setup = "fb_12",
               donor_fe = "idr",
               trial_fe = "idr",
               txd_fe = "idr",
               weight_fe = "sidr",
               weight_is_nested = TRUE,
               censor = 0.8),
    
    data.table(description = "FB Trial 1+2, FEs IDR, Weight nested", # 8
               patch_name = "scen-8-1",
               setup = "fb_12",
               donor_fe = "idr",
               trial_fe = "idr",
               txd_fe = "idr",
               weight_fe = "sidr",
               weight_is_nested = TRUE,
               censor = 0.6),
    
    fill = TRUE
)



# Common options
protocol[, `:=`(patch_dataset = "fb-fes3",
                traits_source = "posterior"
                patch_type = "mean",
                use_traits = "sit",
                fix_donors = "no_Tsym_survivors",
                use_weight = "log",
                use_grm = "",
                group_effect = 0.1,
                prior__latent_period__type = "Fixed",
                prior__latent_period__val2 = 10,
                expand_priors = 3,
                RP_dist = "gamma",
                sample_states = 100,
                ie_output = "true",
                time_step = 1)]

# expand_priors must be moved to before any other priors
if ("expand_priors" %in% names(protocol)) {
    setcolorder(protocol, "expand_priors", before = str_which(names(protocol), "prior__")[[1]])
}

# Labels
protocol[, label := str_c("s", 1:.N)]


# Append "coverage" or "convergence" to description
protocol[, description := str_c(description, ", ", goal)]


## Add replicates ----
n_replicates <- if (goal == "convergence") 1L else 20L
protocol[, scenario := .I]
protocol <- protocol[rep(1:.N, each = n_replicates)]
protocol[, replicate := 1:.N, scenario]
protocol[, dataset := dataset]


# Fixed parameters ----

# Save params along with protocol so we know the defaults for all entries
params <- make_parameters(model_type = "SEIDR",
                          setup = "fb_12",
                          use_traits = "sit",
                          vars = 1.0,
                          cors = 0,
                          group_layout = "fishboost",
                          group_effect = -1,
                          sim_new_data = "r")

params$trial_fe <- ""
params$donor_fe <- ""
params$txd_fe <- ""
params$group_effect <- -1
params$link_traits <- "sildt"
params$sim_link_trial <- "sildt"
params$sim_link_donor <- "sildt"
params$sim_link_shapes <- "ldt"
params$link_trial <- "sildt"
params$link_donor <- "sildt"
params$link_shapes <- "ldt"
params$pass_events <- "Tsym,Tdeath"
params$seed <- if (goal == "convergence") 0 else -1
params$nchains <- if (goal == "convergence") 16 else 4
params$phi <- 1.0
nsample <- 5e6L
params$nsample <- as.integer(nsample)
params$burnin <- as.integer(nsample / 5)
params$thin <- as.integer(max(nsample / 1e4, 1))
params$nsample_per_gen <- as.integer(max(nsample * 2.5 / 1e3))
params$anneal <- "on"
params$anneal_power <- 4
params$sire_version <- "sire22"

# Add missing columns ----

replace_NAs <- function(col) {
    # Make sure columns exist and overwrite NA values with something useful
    if (col %in% names(protocol)) {
        protocol[is.na(get(col)), (col) := params[[col]]]
    } else {
        protocol[, (col) := params[[col]]]
    }
}

missing_cols <- c("model_type", "dataset", "name", "use_traits", "vars", "cors",
                  "sim_new_data", "setup", "group_effect",
                  "trial_fe", "donor_fe", "txd_fe", "weight_is_nested",
                  "sire_version", "nsample_per_gen", # "anneal", "anneal_power", "seed",
                  "phi", "nsample", "burnin", "thin", "nchains")

walk(missing_cols, replace_NAs)

# Tidy up ----

# Give default values for data name, setup, cors, group_layout, ...
protocol[, name := str_c("scen-", scenario, "-", replicate)]

# Prefer to have these columns in this order at the start
cols <-  c("dataset", "description", "scenario", "replicate", "label",
           "name", "sim_new_data", "model_type", "setup")
setcolorder(protocol, intersect(cols, names(protocol)))

message(str_glue("protocol file '{dataset}' has {nrow(protocol)} rows x {ncol(protocol)} cols"))


# Save to file ----
saveRDS(list(protocol = protocol,
             params = params),
        file = str_glue("param_sets/{dataset}.rds"))

# fwrite(protocol, file = "protocol-sim.tsv", sep = "\t", quote = TRUE)
