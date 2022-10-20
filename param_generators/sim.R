library(data.table) # used for all data table manipulations
library(glue)

source("make_parameters.R")

data_set <- "sim-Gsir_cov0_Dlidr-2"

# defaults
{
    default_model <- "SEIDR"
    default_setup <- "fishboost"
    default_use_traits <- "none"
    default_vars <- 1.0
    default_covars <- "none"
    default_group_effect <- -1
    default_donor_fe <- "none"
    default_trial_fe <- "none"
    default_group_layout <- "fishboost"
}


protocol <- rbind(
    data.table(description = "Sim FB Trial 1, SEIDR, traits SIR, cov 0, donor LIDR, large FEs, fixed LP",
               setup = "fb1",
               label = "s1a",
               use_traits = "sir",
               link_traits = "slidr",
               covars = "none",
               donor_fe = "lidr",
               sim_link_donor = "slidr",
               link_donor = "slidr",
               fe_donor_latency = -3,
               fe_donor_infectivity = 2,
               fe_donor_detectability = 3,
               fe_donor_recoverability = 2,
               latent_period = 3,
               prior__latent_period__1 = 2.99,
               prior__latent_period__2 = 3.01,
               group_effect = -1),
    
    data.table(description = "Sim FB Trial 1, SEIDR, traits I, cov 0, donor LIDR, large FEs, fixed LP",
               setup = "fb1",
               label = "s1b",
               use_traits = "i",
               link_traits = "slidr",
               covars = "none",
               donor_fe = "lidr",
               sim_link_donor = "slidr",
               link_donor = "slidr",
               fe_donor_latency = -3,
               fe_donor_infectivity = 2,
               fe_donor_detectability = 3,
               fe_donor_recoverability = 2,
               latent_period = 3,
               prior__latent_period__1 = 2.99,
               prior__latent_period__2 = 3.01,
               group_effect = -1),
    
    fill = TRUE
)

# Save params along with protocol so we know what the defaults were
params <- make_parameters(model_type = "SEIDR",
                          setup = default_setup,
                          use_traits = default_use_traits,
                          vars = default_vars,
                          covars = default_covars,
                          group_layout = default_group_layout,
                          group_effect = default_group_effect,
                          use_fb_data = FALSE)

# cols <- names(protocol)[protocol[1, is.na(.SD)]]
# if (length(cols) > 0) {
#     cols[cols == "use_traits"] <- "traits"
#     protocol[1, cols] <- params[cols]
# }


# Add replicates
n_replicates <- 20L
protocol[, scenario := .I]
protocol <- protocol[rep(1:.N, each = n_replicates)]
protocol[, replicate := 1:n_replicates, by = scenario]
protocol[, data_set := data_set]

# Give default values for data name, setup, covars, group_layout, ...
protocol[, name := paste0("scen-", scenario, "-", replicate)]


# These are things we really want to see in the protocol file, even if they're
# not defined above, probably because the don't have good defaults

replace_NAs <- function(col, val = "") {
    # Make sure columns exist and overwrite NA values with something useful
    if (col %in% names(protocol)) {
        protocol[is.na(get(col)), (col) := val]
    } else {
        protocol[, (col) := val]
    }
}

replace_NAs("model", "SEIDR")
replace_NAs("data_set", data_set)
replace_NAs("name", name)
replace_NAs("use_traits", default_use_traits)
replace_NAs("vars", default_vars)
replace_NAs("covars", default_covars)
replace_NAs("setup", default_setup)
replace_NAs("group_layout", default_group_layout)
replace_NAs("use_fb_data", FALSE)
replace_NAs("trial_fe", "none")
replace_NAs("donor_fe", "none")
# replace_NAs("group_effect", -1)
# replace_NAs("link_traits", "srirr")
# replace_NAs("sim_link_trial", "slidr")
# replace_NAs("sim_link_donor", "slidr")
# replace_NAs("sim_link_shapes", "ldr")
# replace_NAs("link_trial", "slidr")
# replace_NAs("link_donor", "slidr")
# replace_NAs("link_shapes", "ldr")
# replace_NAs("pass_events", 2)
# replace_NAs("seed", 1)
replace_NAs("nsample", 1e7L)
replace_NAs("burnin", 2.5e6L)
replace_NAs("thin", 1e3L)
replace_NAs("quench", "on")
replace_NAs("quench_power", 4)

# Prefer to have these columns in this order at the start
cols <-  c("data_set", "description", "scenario", "replicate", "label",
           "name", "use_fb_data", "model", "setup")
# "trial_fe", "donor_fe",
# "use_traits", "group_effect", "vars", "covars", "group_layout",
# "pass_events", "link_traits", "seed",
# "eta_shape", "rho_shape", "gamma_shape",
# "nsample", "burnin", "thin", "quench", "quench_power")
setcolorder(protocol, intersect(cols, names(protocol)))

message(glue("protocol file '{data_set}' has {nrow(protocol)} rows x {ncol(protocol)} cols"))

protocol_list <- list(protocol = protocol,
                      params = params)

saveRDS(protocol_list, file = glue("data_sets/{data_set}.rds"))

fwrite(protocol, file = "protocol-sim.tsv", sep = "\t", quote = TRUE)
