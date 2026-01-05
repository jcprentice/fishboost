# Common parameters forming a typical base parameter set
# common2 <- dplyr::lst(
common2 <- list(
    model_type = "SEIDR", sim_new_data = "no",
    setup = "fb_12_rpw", use_grm = "pedigree", group_layout = "fishboost",
    link_traits = "sildt", sim_link_traits = "sildt",
    vars = 1, cors = 0.2, group_effect = -1,
    link_shapes = "ldt", sim_link_shapes = "ldt",
    trial_fe  = "",      link_trial  = "sildt", sim_link_trial  = "sildt",
    donor_fe  = "",      link_donor  = "sildt", sim_link_donor  = "sildt",
    txd_fe    = "",      link_txd    = "sildt", sim_link_txd    = "sildt",
    weight_fe = "sildt", link_weight = "sildt", sim_link_weight = "sildt",
    use_weight = "log", weight_is_nested = TRUE,
    cov_prior = "jeffreys", single_prior = "inverse",
    pass_events = "Tsym,Tdeath", time_step = 1,
    # seed = if (goal == "convergence") 0 else -1,
    nchains = if (goal == "convergence") 16 else 4,
    nsample = 1e6, burnprop = 0.2, thinto = 1e4,
    nsample_per_gen = 1, phi = 1.0,
    anneal = "on", anneal_power = 4, sire_version = "bici", bici_cmd = "inf"
)

common2$nsample_per_gen <- max(3e-3 * common2$nsample, 1)

# Order of first few columns
cols <-  c("dataset", "description", "scenario", "replicate", "label",
           "name", "sim_new_data", "model_type", "setup")

