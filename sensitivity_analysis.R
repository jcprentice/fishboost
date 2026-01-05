# Load Libraries ----
{
    library(data.table) # used for all data table manipulations
    library(purrr)      # functional programming
    library(Matrix)     # allows Matrix
    library(xml2)       # used to generate XML file
    library(corpcor)    # efficient estimation of Cov and (partial) Corr
    library(MASS)       # needed for mvnorm
    library(MCMCglmm)
    library(ggplot2)    # plotting
    library(tictoc)     # check performance
    library(codetools)  # check for global variables
    library(here)       # useful according to Hadley Wickham
}

# Load source files ----
{
    source("make_parameters.R")
    source("get_R0.R")
    source("make_pedigree.R")
    source("make_grm.R")
    source("make_traits_from_grm.R")
    source("make_traits_from_pedigree.R")
    source("set_groups.R")
    source("simulate_epidemic.R")
    source("init_populations.R")
    source("prepare_data.R")
    source("generate_sire_xml.R")
    source("make_time_series.R")
    source("make_plots.R")
    source("utils.R")
}

# Helper functions ----
epi_peak_I <- function(X) X[, max(I)]
epi_peak_t <- function(X) X[, time[which.max(I)]]
epi_duration <- function(X) X[.N, time]

# SA for eta ----
{
    # Set up initial population structure ----
    params_seir <- make_parameters(model_type = "SEIR", setup = "groups")
    params_sir  <- make_parameters(model_type = "SIR", setup = "groups")
    pedigree    <- make_pedigree(params_seir)
    
    
    # List of parameters to simulate over ----
    num_pts <- 51L
    latent_periods <- seq(0, 50, length.out = num_pts)
    
    # Initialise results DT
    results1 <- data.table(
        id = seq_len(num_pts),
        model = factor(c("sir", rep("seir", num_pts - 1))),
        latency = latent_periods,
        R0 = 0.0,
        peak_I = 0L,
        peak_t = 0.0,
        duration = 0.0
    )
    
    
    for (i in 1:num_pts) {
        message(str_glue("\nIteration {i} of {num_pts}"))
        
        params <- if (i == 1) params_sir else params_seir
        
        params$LP_scale <- latent_periods[i]
        popn <- pedigree |>
            make_traits_from_pedigree(params) |>
            set_groups(params) |>
            simulate_epidemic(params)
        
        ts <- make_time_series(popn, params)
        
        results1[id == i,
                 `:=`(R0 = get_R0(popn),
                      peak = epi_peak(ts),
                      duration = epi_duration(ts))]
        set(results1, i, c("R0", "peak_I", "peak_t", "duration"),
            list(get_R0(popn), epi_peak_I(ts), epi_peak_t(ts), epi_duration(ts)))
    }
    
    print(results1)
    
    # Plot the results ----
    plt1 <- ggplot(results1, aes(x = 1 / eta)) +
        scale_y_continuous(trans = "log10") +
        geom_point(aes(y = peak_I, colour = "peak infections")) +
        geom_point(aes(y = peak_t, colour = "time of peak")) +
        geom_point(aes(y = R0, colour = "R0")) +
        geom_point(aes(y = duration, colour = "duration")) +
        geom_smooth(aes(y = peak_I, colour = "peak infections")) +
        geom_smooth(aes(y = peak_t, colour = "time of peak")) +
        geom_smooth(aes(y = R0, colour = "R0")) +
        geom_smooth(aes(y = duration, colour = "duration")) +
        labs(x = "Latent period (days)", y = NULL)
    
    print(plt1)
    
    saveRDS(mget(c("params_seir", "params_sir", "popn", "results1", "plt1")),
            file = "data/sa1.rds")
}


# Test if SEIR = SIR when latent period = 0 ----
{
    # set up initial populations
    params_seir <- make_parameters(model_type = "SEIR", setup = "")
    params_seir$LP_scale <- 0
    params_sir  <- make_parameters(model_type = "SIR", setup = "")
    popn <- make_pedigree(params_seir) |>
        make_traits_from_pedigree(params_seir) |>
        set_groups(params_seir)
    
    # storage for results
    models <- c("sir", "seir")
    n <- 1000
    results2 <- data.table(expand.grid(id = 1:n, model = models))
    
    for (model_type in models) {
        params <- get(str_c("params_", model_type))
        for (i in 1:n) {
            message(str_glue("Itn {i} / {n}"))
            popn2 <- simulate_epidemic(popn, params)
            ts <- make_time_series(popn2, params)
            
            results2[id == i & model == model_type,
                     `:=`(R0 = get_R0(popn2),
                          peak_I = epi_peak_I(ts),
                          peak_t = epi_peak_t(ts),
                          duration = epi_duration(ts))]
        }
    }
    
    means <- results2[, map(.SD, mean), by = model, .SDcols = -"id"]
    print(means)
    
    plt2 <- ggplot(results2, aes(x = peak_t, fill = model)) +
        geom_histogram(alpha = 0.5, position = "identity", bins = 30) +
        labs(x = "Time of peak (days)")
    
    print(plt2)
    
    # Save SA data ----
    # This data takes a while to generate and is worth saving
    saveRDS(mget(c("params_seir", "params_sir", "popn", "results2", "plt2")),
            file = "data/sa2.rds")
}

# Sensitivity Analysis for covariance
{
    params   <- make_parameters(model_type = "SIR", setup = "chris")
    pedigree <- make_pedigree(params)
    
    Sigma_G <- params$Sigma_G
    Sigma_E <- params$Sigma_E
    
    num_pts <- 51
    cov_mult <- seq(0, 2, length.out = num_pts + 1)[-1]
    
    results3 <- data.table(id = 1:num_pts, cov_mult = cov_mult)
    
    for (i in 1:num_pts) {
        params$Sigma_E <- Sigma_E * cov_mult[i]
        params$Sigma_G <- Sigma_G * cov_mult[i]
        
        popn <- make_traits_from_pedigree(params) |>
            set_groups(params) |>
            simulate_epidemic(popn, params)
        
        ts <- make_time_series(popn, params)
        
        results3[id == i, `:=`(R0 = get_R0(popn),
                               peak_I = epi_peak_I(ts),
                               peak_t = epi_peak_t(ts),
                               duration = epi_duration(ts))]
    }
    
    # Plot the results ----
    plt3 <- ggplot(results3, aes(x = cov_mult)) +
        scale_y_continuous(trans = "log10") +
        geom_point(aes(y = peak_I, colour = "peak infections")) +
        geom_point(aes(y = peak_t, colour = "time of peak")) +
        geom_point(aes(y = R0, colour = "R0")) +
        geom_point(aes(y = duration, colour = "duration")) +
        geom_smooth(aes(y = peak_I, colour = "peak infections")) +
        geom_smooth(aes(y = peak_t, colour = "time of peak")) +
        geom_smooth(aes(y = R0, colour = "R0")) +
        geom_smooth(aes(y = duration, colour = "duration")) +
        labs(x = "Multiplier for covariance", y = NULL)
    
    print(plt3)
    
    saveRDS(mget(c("params_sir", "popn", "results3", "plt3")),
            file = "data/sa3.rds")
}
