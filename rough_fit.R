source("libraries.R")
source("source_files.R")

library(cowplot)

load("data/fb_plots.RData")

params <- make_parameters(
    model_type = "SEIDR", # "SIR", "SEIR", "SIDR", or "SEIDR"
    name = "expt1",
    setup = "fb1", # chris, small, fishboost, fb1, fb2
    use_traits = "sir", # "all", "none", "sir", "si" etc.
    vars = 0.5, # c(0.5, 1.2, 0.5)
    covars = "none", # 0.2, "positive", "negative", "mixed", "sir_only" etc.
    group_layout = "fishboost", # "random", "family", "striped", "fishboost"
    trial_fe = "",
    donor_fe = "l",
    txd_fe = "",
    use_fb_data = FALSE # don't simulate, just use FB data
)

params$show_plots <- FALSE
params$tmax <- 102 # 159

x <- with(params,
          round(c(beta=r_beta,
                  LP=1/r_eta, eta_shape=r_eta_shape,
                  DP=1/r_rho, rho_shape=r_rho_shape,
                  RP=1/r_gamma, gamma_shape=r_gamma_shape),
                digits=2))
print(x)

params$fe_vals["donor", "latency"] <- -1
params$fe_vals["donor", "detectability"] <- -1

{
    params$r_beta <- 0.08
    LP <- 3.14 / exp(0.8 * params$fe_vals["donor", "latency"])
    params$r_eta <- 1 / LP
    params$r_eta_shape <- 2
    params$r_eta_rate <- params$r_eta_shape / LP
    DP <- 3.14 / exp(0.8 * params$fe_vals["donor", "detectability"])
    params$r_rho <- 1 / DP
    params$r_rho_shape <- 2
    params$r_rho_rate <- params$r_rho_shape / DP
    RP <- 8.19
    params$r_gamma <- 1 / RP
    params$r_gamma_shape <- 1.22
    params$r_gamma_rate <- params$r_rho_shape / RP
}

{
    params$group_effect <- -1
    # Quick check of how things are looking
    summarise_params(params)

    pedigree <- make_pedigree(params)
    GRM <- make_grm(pedigree)
    traits <- make_traits_from_pedigree(pedigree, params)
    traits <- set_groups(traits, params)
    traits <- apply_fixed_effects(traits, params)

    pop <- simulate_epidemic(traits, params)
    params$estimated_R0 <- get_R0(pop)

    plt <- plot_SxxDR(pop, params)

    events <- make_time_series_seidr(pop, params)
    t_end <- events[time <= params$tmax][time == max(time)]
    print(t_end)

    gc()
    plt2 <- plot_grid(plotlist = list(plt_fb1, plt))
    print(plt2)
}

#     time   S  E I  D   R ID
# FB1: 102 175  1 0 39 685 39
# FB2: 159 250 53 0 16 581 16

#    parameter      mean
# ----------------------
#         beta  0.193852
#           LP  2.696651
#    eta shape  1.000225
#           DP  6.377101
#    rho shape  0.415509
#           RP 11.973963
#  gamma shape  0.831726
#        sigma  0.119679

