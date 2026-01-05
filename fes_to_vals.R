{
    library(purrr)
    library(data.table)
    library(stringr)
    source("utils.R")
}

fes_to_vals <- function(dataset = "fb-lp", scen = 1, rep = 1) {
    # dataset <- "fb-lp"; scen <- 1; rep <- 1
    f <- str_glue("datasets/{dataset}/results/scen-{scen}-{rep}.rds")
    if (!file.exists(f)) return(NULL)
    res <- readRDS(f)
    params <- res$params
    pe <- res$parameter_estimates
    fe <- pe[str_detect(parameter, "trial|donor|txd"),
             setNames(mean, parameter)] |>
        as.list()
    
    fes <- expand.grid(c("trial", "donor", "txd"),
                       c("s", "l", "i", "d", "r")) |>
        apply(1, str_flatten, "_")
    
    walk(fes, \(x) fe[[x]] <<- fe[[x]] %||% 0)
    
    pop <- res$popn[sdp == "progeny", .N, .(trial, donor)]
    pop[, `:=`(group = str_c("Tr", trial, ",", fifelse(donor == 1, "Don", "Rec")),
               p = N / sum(N),
               N = NULL)]
    
    pop[, `:=`(trial_p = (trial == 2) - sum(p * (trial == 2)),
               donor_p = (donor == 1) - sum(p * (donor == 1)),
               txd_p = (donor == 1 & trial == 2) - sum(p * (donor == 1 & trial == 2)),
               beta = params$r_beta,
               LP = params$latent_period,
               DP = params$detection_period,
               RP = params$removal_period)]
    
    pop[, `:=`(
        # beta_ = beta * exp(trial_p * fe$trial_i + donor_p * fe$donor_i + txd_p * fe$txd_i),
        beta_ = beta * exp(trial_p * fe$trial_i),
        LP_   = LP   * exp(trial_p * fe$trial_l + donor_p * fe$donor_l + txd_p * fe$txd_l),
        DP_   = DP   * exp(trial_p * fe$trial_d + donor_p * fe$donor_d + txd_p * fe$txd_d),
        RP_   = RP   * exp(trial_p * fe$trial_t + donor_p * fe$donor_t + txd_p * fe$txd_t)
    )]
    
    vals <- pop[, .(group, beta = round(beta_, 2),
                    latent_period = round(LP_, 1),
                    detection_period = round(DP_, 1),
                    removal_period = round(RP_, 1))] |>
        melt(id.var = "group")
    vals <- vals[-c(2, 4), .(parameter = str_c(variable, "_", group), true_val = value)]
    
    vals[str_detect(parameter, "beta"),
         parameter := str_remove(parameter, ",Don")]
}

# vals_to_fes <- function(dataset = "fb-test", scen = 1) {
#     # dataset <- "fb-test"; scen <- 1
#     f <- str_glue("datasets/{dataset}/results/scen-{scen}-1.rds")
#     params <- readRDS(f)$params
#     pop <- readRDS(f)$popn[sdp == "progeny", .N, .(trial, donor)]
#     pe <- readRDS(f)$parameter_estimates
#     fe <- pe[str_detect(parameter, "beta|period"), setNames(mean, parameter)] |>
#         as.list()
#     
#     pop[, `:=`(group = str_c("Tr", trial, ",", fifelse(donor == 1, "Don", "Rec")),
#                p = N / sum(N),
#                N = NULL)]
#     
#     pop[, `:=`(trial_p = (trial == 2) - sum(p * (trial == 2)),
#                donor_p = (donor == 1) - sum(p * (donor == 1)),
#                txd_p = (donor == 1 & trial == 2) - sum(p * (donor == 1 & trial == 2)),
#                beta = params$r_beta,
#                LP = params$latent_period,
#                DP = params$detection_period,
#                RP = params$removal_period)]
# }
