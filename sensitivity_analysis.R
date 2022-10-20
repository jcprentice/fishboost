# Load Libraries ----
{
    library(data.table) # used for all data table manipulations
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
    # source("discretise.R")
    source("prepare_data.R")
    source("generate_sire20_xml.R")
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
    latency_periods <- seq(0, 50, length.out = num_pts)
    etas <- 1 / latency_periods

    # Initialise results DT
    results1 <- data.table(
        id = seq(num_pts),
        model = factor(c("sir", rep("seir", num_pts - 1))),
        latency = latency_periods,
        eta = etas,
        R0 = 0.0,
        peak_I = 0L,
        peak_t = 0.0,
        duration = 0.0
    )


    for (i in 1:num_pts) {
        message("\nIteration ", i, " of ", num_pts)

        params <- if (i == 1) params_sir else params_seir

        params$r_eta_rate <- etas[i]
        traits  <- make_traits_from_pedigree(pedigree, params)
        traits1 <- set_groups(traits, params)

        pop     <- simulate_epidemic(traits1, params)

        ts <- make_time_series(pop, params)

        results1[id == i,
                `:=`(R0 = get_R0(pop),
                     peak = epi_peak(ts),
                     duration = epi_duration(ts))]
        set(results1, i, c("R0", "peak_I", "peak_t", "duration"),
            list(get_R0(pop), epi_peak_I(ts), epi_peak_t(ts), epi_duration(ts)))
    }

    print(results1)

    # Plot the results ----
    plt1 <- ggplot(results1, aes(x = 1/eta)) +
        scale_y_continuous(trans = "log10") +
        geom_point(aes(y = peak_I, colour = "peak infections")) +
        geom_point(aes(y = peak_t, colour = "time of peak")) +
        geom_point(aes(y = R0, colour = "R0")) +
        geom_point(aes(y = duration, colour = "duration")) +
        geom_smooth(aes(y = peak_I, colour = "peak infections")) +
        geom_smooth(aes(y = peak_t, colour = "time of peak")) +
        geom_smooth(aes(y = R0, colour = "R0")) +
        geom_smooth(aes(y = duration, colour = "duration")) +
        labs(x = "Latency period (days)", y = NULL)

    print(plt1)

    save(params_seir, params_sir, traits1, results1, plt1, file = "data/sa1.RData")
}


# Test if SEIR = SIR when latency period = 0 ----
{
    # set up initial populations
    params_seir <- make_parameters(model_type = "SEIR", setup = "")
    params_seir$r_eta_rate <- Inf
    params_sir  <- make_parameters(model_type = "SIR", setup = "")
    pedigree <- make_pedigree(params_seir)
    traits <- make_traits_from_pedigree(pedigree, params_seir)
    traits2 <- set_groups(traits, params_seir)

    # storage for results
    models <- c("sir", "seir")
    n <- 1000
    results2 <- data.table(expand.grid(id = 1:n, model = models))

    for (model_type in models) {
        params <- get(paste0("params_", model_type))
        for (i in 1:n) {
            message("Itn ", i, " / ", n)
            pop <- simulate_epidemic(traits2, params)
            ts <- make_time_series(pop, params)

            results2[id == i & model == model_type,
              `:=`(R0 = get_R0(pop),
                   peak_I = epi_peak_I(ts),
                   peak_t = epi_peak_t(ts),
                   duration = epi_duration(ts))]
        }
    }

    means <- results2[, lapply(.SD, mean), by = model, .SDcols = -"id"]
    print(means)

    plt2 <- ggplot(results2, aes(x = peak_t, fill = model)) +
        geom_histogram(alpha = 0.5, position = "identity", bins = 30) +
        labs(x = "Time of peak (days)")

    print(plt2)

    # Save SA data ----
    # This data takes a while to generate and is worth saving
    save(params_seir, params_sir, traits2, results2, plt2, file = "data/sa2.RData")
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

        traits <- make_traits_from_pedigree(pedigree, params)
        traits3 <- set_groups(traits, params)

        pop <- simulate_epidemic(traits3, params)
        ts <- make_time_series(pop, params)

        results3[id == i, `:=`(R0 = get_R0(pop),
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

    save(params_sir, traits3, results3, plt3, file = "data/sa3.RData")
}
