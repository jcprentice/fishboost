library(data.table)

# Generate parameter list
make_parameters <- function(
    # epidemic model: "SIR", "SIDR", "SEIR", "SEIDR", "...-P" for environmental
    model_type = "SIDR",
    # filename for data
    name = "foo",
    # folder suffix for data (default "")
    data_set = "",
    # population layout
    setup = "fishboost",
    # which traits to use (this is going to be "clever")
    # use "all" or "none", or a subset using the first letter of each trait
    # "slidr" = "sus", "lat", "inf", "det", "rec"
    use_traits = "sir",
    # variances: diagonal of cov matrix, either a single value, or a vector, in
    # which case it should be compatible with "use_traits"
    vars = 0.5,
    # cov matrix, numerical: 0.2, c(0.2, 0.2), or character: "chris", "positive", "sir_only"
    covars = 0.2,
    # how to arrange individuals into groups: "random", "striped", "family", "fishboost"
    group_layout = "random",
    # simulate new data or use Fishboost data
    use_fb_data = FALSE
) {

    message("Setting parameters ...")

    # Fix input NAs ----

    # check if any values are NA or NULL, and replace with sensible defaults
    get_default <- function(x, x_def) if (is.na(x) || is.null(x)) x_def else x

    # uncomment if we need these variables
    # model_type <- "SIR"; name <- "foo"; data_set <- ""; setup <- "fishboost"; use_traits <- "all"; vars <- 0.5; covars <- 0.2; group_layout <- "random"; use_fb_data <- FALSE;

    model_type   <- get_default(model_type, "SIR")
    name         <- get_default(name, paste0(tolower(model_type), 1L))
    data_set     <- get_default(data_set, "")
    setup        <- get_default(setup, "")
    use_traits   <- get_default(use_traits, "all")
    vars         <- get_default(vars, 0.5)
    covars       <- get_default(covars, 0.2)
    group_layout <- get_default(group_layout, "random")
    use_fb_data  <- get_default(use_fb_data, FALSE)


    # Set output directories ----

    # Directories should be something like "data", "data-sim-xyz", or
    # "data-fb-xyz", so check if data_set starts with "-" (unless it's empty).
    if (grepl("^-", data_set) == FALSE && data_set != "") {
        data_set <- paste0("-", data_set)
    }
    data_dir    <- glue("data{data_set}")
    results_dir <- sub("data", "results", data_dir)


    # Population setup ----
    message(glue(" - Setup is: '{setup}'\n",
                 " - Group layout is: '{group_layout}'"))

    # note: in balanced populations, dpsire = dams per sire, ppdam = progeny per dam
    switch(setup,
           "fishboost" = {
               nsires <- 29L; ndams <- 25L; nprogeny <- 1800L;
               dpsire <- 1L; ppdam <- 72L;
               ngroups <- 72L; I0 <- 5L;
           }, "fb1" = {
               nsires <- 14L; ndams <- 14L; nprogeny <- 900L;
               dpsire <- 1L; ppdam <- 72L;
               ngroups <- 36L; I0 <- 5L;
           }, "fb2" = {
               nsires <- 18L; ndams <- 14L; nprogeny <- 900L;
               dpsire <- 1L; ppdam <- 72L;
               ngroups <- 36L; I0 <- 5L;
           }, "chris" = {
               nsires <- 100L; ndams <- 2000L; nprogeny <- 2000L;
               dpsire <- 20L; ppdam <- 1L;
               ngroups <- 200L; I0 <- 1L;
           }, "small" = {
               nsires <- 3L; ndams <- 6L; nprogeny <- 12L;
               dpsire <- 2L; ppdam <- 2L;
               ngroups <- 4L; I0 <- 1L;
           }, {
               nsires <- 3L; ndams <- 6L; nprogeny <- 12L;
               dpsire <- 2L; ppdam <- 2L;
               ngroups <- 4L; I0 <- 1L;
           }
    )

    # derived numbers
    nparents <- nsires + ndams
    ntotal <- nprogeny + nparents
    group_size <- nprogeny / ngroups


    # Traits ----
    message(glue(" - Model type is: '{model_type}'"))

    # traitnames        list of traits that the model uses
    # corr_signs        for when covariance is mixed, should more disease = higher (+1) or lower (-1)
    # compartments      names of compartments
    # timings           name for time when individual enters compartment
    # use_parameters    subset of parameters to pass to SIRE 2.1

    # check for reservoir model
    reservoir <- endsWith(model_type, "-res")
    if (reservoir) {
        model_type <- sub("-res", "", model_type)
    }

    # compartments are just the letters in model_type (ignore repeated S)
    compartments <- unique(unlist(strsplit(model_type, "")))

    switch(
        model_type,
        "SEIDR" = {
            all_traitnames <- c("susceptibility", "latency", "infectivity", "detectability", "recoverability")
            corr_signs <- c(1, -1, 1, 1, -1)
            timings <- c("Tinf", "Tinc", "Tsym", "Trec")
            use_parameters <- c("beta", "latent_period", "eta_shape", "detection_period",
                                "rho_shape", "recovery_period", "gamma_shape")
        }, "SIDR" = {
            all_traitnames <- c("susceptibility", "infectivity", "detectability", "recoverability")
            corr_signs <- c(1, 1, 1, -1)
            timings <- c("Tinf", "Tsym", "Trec")
            use_parameters <- c("beta", "detection_period", "rho_shape", "recovery_period", "gamma_shape")
        }, "SEIR" = {
            all_traitnames <- c("susceptibility", "latency", "infectivity", "recoverability")
            corr_signs <- c(1, -1, 1, -1)
            timings <- c("Tinf", "Tsym", "Trec")
            use_parameters <- c("beta", "latent_period", "eta_shape", "recovery_period", "gamma_shape")
        }, "SIR" = {
            all_traitnames <- c("susceptibility", "infectivity", "recoverability")
            corr_signs <- c(1, 1, -1)
            timings <- c("Tinf", "Trec")
            use_parameters <- c("beta", "recovery_period", "gamma_shape")
        }, "SIS" = {
            all_traitnames <- c("susceptibility", "infectivity", "recoverability")
            corr_signs <- c(1, 1, -1)
            timings <- c("Tinf", "Trec")
            use_parameters <- c("beta", "recovery_period", "gamma_shape")
        }, "SI" = {
            all_traitnames <- c("susceptibility", "infectivity")
            corr_signs <- c(1, 1)
            timings <- c("Tinf")
            use_parameters <- c("beta")
        }, {
            stop(" - Unknown model!")
        }
    )

    ## Genetic and Environmental covariances ----

    # We might be given a subset of traits to use
    traitnames <- switch(
        use_traits,
        "all" = all_traitnames,
        "none" = c(),
        "seir" = c("susceptibility", "latency", "infectivity", "recoverability"),
        "sidr" = c("susceptibility", "infectivity", "detectibility", "recoverability"),
        "sir" = c("susceptibility", "infectivity", "recoverability"),
        "si" = c("susceptibility", "infectivity"),
        stop(" - Unknown 'traitnames'!")
    )
    ntraits <- length(traitnames)

    message(" - ", ntraits, " traits: ", paste0(traitnames, collapse = ", "))

    # Assume default of var(i) = 0.5 and cov(i, j) = 0
    Sigma_G <- diag(vars[vars > 0], ntraits)
    dimnames(Sigma_G) <- list(traitnames, traitnames)
    Sigma_E <- Sigma_G


    # Useful to have these indices (since they depend on the model). I could
    # just write out the trait name each time, but this is easier!
    iS <- which(traitnames == "susceptibility")
    iI <- which(traitnames == "infectivity")
    iR <- which(traitnames == "recoverability")

    sir <- c(iS, iI, iR)
    si  <- c(iS, iI)

    # covars is either a vector (possibly of length 1) or a string
    if (is.numeric(covars)) {
        Sigma_G[upper.tri(Sigma_G)] <- covars
        Sigma_G[lower.tri(Sigma_G)] <- t(Sigma_G)[lower.tri(Sigma_G)]
        Sigma_E <- Sigma_G
    } else {
        switch(
            covars,
            "chris" = {
                diag(Sigma_G) <- 1e-6;
                diag(Sigma_E) <- 1e-6;
                Sigma_G[sir, sir] <- c(+0.33,  +0.3,  -0.044,
                                       +0.3,   +1.68, -0.5,
                                       -0.044, -0.5,  +0.6);
                Sigma_E[sir, sir] <- c(+0.77,  +0.56, -0.25,
                                       +0.56,  +1.12, -0.4,
                                       -0.25,  -0.4,  +0.9)
            }, "positive" = {
                Sigma_G[] <- 0.2
                diag(Sigma_G) <- 0.5
            }, "negative" = {
                Sigma_G[] <- -0.2
                diag(Sigma_G) <- 0.5
            }, "mixed" = {
                Sigma_G[] <- 0.2
                diag(Sigma_G) <- 0.5
                Sigma_G <- Sigma_G * corr_signs %*% t(corr_signs)
            }, "sir_no_cov" = {
                diag(Sigma_G) <- 1e-6
                diag(Sigma_G[sir, sir]) <- 0.5
            }, "sir_pos_cov" = {
                Sigma_G[sir, sir] <- 0.2
                diag(Sigma_G) <- 1e-6
                diag(Sigma_G[sir, sir]) <- 0.5
            }, "si_no_cov" = {
                diag(Sigma_G) <- 1e-6
                diag(Sigma_G[si, si]) <- 0.5
            }, "si_pos_cov" = {
                Sigma_G[si, si] <- 0.2
                diag(Sigma_G) <- 1e-6
                diag(Sigma_G[si, si]) <- 0.5
            }, "sir_only" = {
                diag(Sigma_G) <- 1e-6
                diag(Sigma_G[sir, sir]) <- 0.5
            } # none is default and already done
        )
    }

    message(" - Covariance is: '", paste0(covars, collapse = ", "), "'")

    # Sigma_G and Sigma_E are the same in all cases except for Chris's
    if (covars != "chris") {
        Sigma_E <- Sigma_G
    }


    ### Heritability and correlation ----

    h2 <- Sigma_G / (Sigma_G + Sigma_E)
    if (ntraits > 0) {
        cor_G <- cov2cor(Sigma_G)
        cor_E <- cov2cor(Sigma_E)
    } else {
        cor_G <- matrix(0, 0, 0)
        cor_E <- matrix(0, 0, 0)
    }


    # Epidemic parameters ----

    ## Infection coefficient ----
    r_beta <- 0.1

    ## Reservoir ----
    r_zeta   <- 0.05 # shedding
    r_lambda <- 0.5 # decay

    # Note for Gamma dist:
    # shape = mean^2 / var
    # rate  = mean / var

    ## Latent Periods ----

    latent_period <-  6.282313

    ### Incubation period ----

    # Calculated for an SEIR model based on donors in trial 1 (+ 0.5)
    # fitdist(Tsym + 0.5, "gamma")
    r_eta <- 1.0 / latent_period
    r_eta_shape <- 1.4565712 # 2.0 # 10.0
    r_eta_rate  <- 0.2318177 # r_eta_shape / latent_period


    ### Asymptomatic infectious period ----

    # This is currently identical to the latency
    r_rho <- 1.0 / latent_period
    r_rho_shape <- 1.4565712 # 2.0 # 10.0
    r_rho_rate  <- 0.2318177 # r_eta_shape / latent_period


    # If it's an SEIDR model, we want the combined latent period (S->E->I->D) to
    # have the same distribution for the LP in the SEIR and SIDR models, which
    # we do by halving each shape parameter
    if (model_type == "SEIDR") {
        latent_period <- latent_period / 2
        r_eta <- r_eta / 2
        r_rho <- r_rho / 2
        r_eta_shape <- r_eta_shape / 2
        r_rho_shape <- r_rho_shape / 2
    }


    ## Recovery ----

    # fitdist(RP + 0.5, "gamma")
    recovery_period <- 10.97461
    r_gamma <- 1.0 / recovery_period
    r_gamma_shape <- 0.92780265
    r_gamma_rate <- 0.08455325 # r_gamma_shape / recovery_period


    ## Group effect ----
    group_effect <- -1


    ## R0 ----

    # This should target R0 around 5?
    R0 <- r_beta / (r_gamma + if ("D" %in% compartments) r_rho else 0)

    # SIS models might reach an endemic equilibrium and not terminate, so need
    # to set maximum time to run the simulation (this might need adjusting
    # depending on the other parameters)
    tmax <- if (model_type == "SIS") 20 else NA_real_

    # Time interval if SIRE only sees observations every few days (note: for
    # continuous time, interval = 0)
    t_gap <- 0.0 # 7.0



    # MCMC settings ----

    ## Samples ----
    nsample <- 1e6L
    burnin  <- nsample %/% 5L
    thin    <- max(nsample %/% 1e4L, 1L)

    ## Priors ----

    # Trying to work with the different parameters in SIRE 2.0 and 2.1
    priors <- as.data.table(dplyr::tribble(
        ~parameter, ~symbol, ~val1, ~val2, ~true_val, ~use_s20, ~use,

        "beta",             "β",   0, 1,                    r_beta,          TRUE,  TRUE,
        "latent_period",    "",    0, 20 * latent_period,   latent_period,   FALSE, FALSE,
        "eta_shape",        "",    0, 5 * r_eta_shape,      r_eta_shape,     FALSE, FALSE,
        "detection_period", "",    0, 20 * latent_period,   latent_period,   FALSE, FALSE,
        "rho_shape",        "",    0, 5 * r_rho_shape,      r_rho_shape,     FALSE, FALSE,
        "recovery_period",  "",    0, 20 * recovery_period, recovery_period, FALSE, FALSE,
        "",                 "γ",   0, 2,                    r_gamma,         TRUE,  FALSE,
        "gamma_shape",      "k",   0, 2 * r_gamma_shape,    r_gamma_shape,   TRUE,  FALSE,
        "",                 "G",  -3, 3,                    NA_real_,        TRUE,  FALSE,
        "sigma",            "σ_G", 0, 0.5,                  group_effect,    TRUE,  TRUE,

        "",          "q_g",  -5,    5,   1,               TRUE,  FALSE,
        "",          "q_f",  -5,    5,   1,               TRUE,  FALSE,
        "",          "q_r",  -5,    5,   1,               TRUE,  FALSE,
        "",          "ε_g",  -5,    5,   1,               TRUE,  FALSE,
        "",          "ε_f",  -5,    5,   1,               TRUE,  FALSE,
        "",          "ε_r",  -5,    5,   1,               TRUE,  FALSE,

        "cov_G_ss",  "Ω_gg",  1e-3, 3,    Sigma_G[iS, iS], TRUE,  TRUE,

        "cov_G_ii",  "Ω_ff",  1e-3, 3,    Sigma_G[iI, iI], TRUE,  TRUE,
        "cov_G_rr",  "Ω_rr",  1e-3, 3,    Sigma_G[iR, iR], TRUE,  TRUE,
        "",          "Ω_gf", -1,    1,    Sigma_G[iS, iI], TRUE,  FALSE,
        "",          "Ω_gr", -1,    1,    Sigma_G[iS, iR], TRUE,  FALSE,
        "",          "Ω_fr", -1,    1,    Sigma_G[iI, iR], TRUE,  FALSE,
        "r_G_si",    "",     -0.95, 0.95, cor_G[iS, iI],   FALSE, TRUE,
        "r_G_sr",    "",     -0.95, 0.95, cor_G[iS, iR],   FALSE, TRUE,
        "r_G_ir",    "",     -0.95, 0.95, cor_G[iI, iR],   FALSE, TRUE,

        "cov_E_ss",  "Ψ_gg",  1e-3, 3,    Sigma_E[iS, iS], TRUE,  TRUE,
        "cov_E_ii",  "Ψ_ff",  1e-3, 3,    Sigma_E[iI, iI], TRUE,  TRUE,
        "cov_E_rr",  "Ψ_rr",  1e-3, 3,    Sigma_E[iR, iR], TRUE,  TRUE,
        "",          "Ψ_gf", -1,    1,    Sigma_E[iS, iI], TRUE,  FALSE,
        "",          "Ψ_gr", -1,    1,    Sigma_E[iS, iR], TRUE,  FALSE,
        "",          "Ψ_fr", -1,    1,    Sigma_E[iI, iR], TRUE,  FALSE,
        "r_E_si",    "",     -0.95,  0.95, cor_E[iS, iI],  FALSE, TRUE,
        "r_E_sr",    "",     -0.95,  0.95, cor_E[iS, iR],  FALSE, TRUE,
        "r_E_ir",    "",     -0.95,  0.95, cor_E[iI, iR],  FALSE, TRUE
    ))

    # Some priors set to not use with SIRE 2.1, update if we really do want them
    priors[parameter %in% use_parameters, use := TRUE]
    if (group_effect < 0) {
        priors[parameter == "sigma", use := FALSE]
    }

    # remove priors if traits not used
    if (!1 %in% sir) {
        priors[parameter %in% c("cov_G_ss", "r_G_si", "r_G_sr", "cov_E_ss",
                                "r_E_si", "r_E_sr"), use := FALSE]
    }
    if (!2 %in% sir) {
        priors[parameter %in% c("cov_G_ii", "r_G_si", "r_G_ir", "cov_E_ii",
                                "r_E_si", "r_E_ir"), use := FALSE]
    }
    if (!3 %in% sir) {
        priors[parameter %in% c("cov_G_rr", "r_G_sr", "r_G_ir", "cov_E_rr",
                                "r_E_sr", "r_E_ir"), use := FALSE]
    }

    # Additional settings ----

    # choose between SIRE 2.0 and 2.1
    use_sire_21 <- TRUE

    # What we pass to SIRE 2.0, choice of
    # "Tinf": actual infection time
    # "Tsym": time of symptoms
    # "estimated_Tinf_from_donors": Tsym - mean incubation period of donors
    # "estimated_Tinf_per_individual": Tsym with all gaps covered
    pass_Tsym <- "Tsym"

    # Do we clear Tsym times so SIRE has to guess them?
    clear_Tsym <- FALSE

    # Show details
    DEBUG <- FALSE


    # Create params list ----
    #
    # This will be a list for passing to functions
    params <- mget(
        # inputs
        c("data_set", "name", "data_dir", "results_dir",
          "model_type", "use_fb_data", "setup", "covars", "group_layout",
          # population
          "nsires", "ndams", "nparents", "nprogeny", "ngroups", "ntotal",
          "dpsire", "ppdam", "group_size", "I0",
          # model traits
          "compartments", "reservoir", "timings", "all_traitnames", "traitnames",
          "ntraits", "Sigma_E", "Sigma_G", "cor_G", "cor_E", "h2",
          # model parameters (rates)
          "r_beta", "r_zeta", "r_lambda",                                # infection
          "latent_period", "r_eta", "r_eta_shape", "r_eta_rate",         # latency
          "r_rho", "r_rho_shape", "r_rho_rate",                          # detection
          "recovery_period", "r_gamma", "r_gamma_shape", "r_gamma_rate", # recovery
          "R0", "group_effect",
          # mcmc & extra
          "t_end", "tmax", "t_gap",
          "nsample", "burnin", "thin", "priors",
          "use_sire_21", "pass_Tsym", "clear_Tsym", "DEBUG"))

    params
}


# Handy summary ----
summarise_params <- function(params) {
    with(params, {
        if (use_fb_data) {
            message("Using Fishboost data")
        } else {
            message(glue("Simulating new data with {model_type} model"))
        }

        message(glue(
            " - Demography is:\n",
            "     {nsires} sires, {ndams} dams, {nprogeny} progeny ",
            "({ntotal} total), {ngroups} groups ({group_size} per group)\n",
            "     R0 = {R0}\n)",
            " - Running MCMC with:\n",
            "     It = '{pass_Tsym}'\n",
            "     {nsample} samples / {burnin}, {burnin} / {thin} thin\n",
            " - Data will be saved to '{data_dir}/{name}.xml'"
        ))
    })
}
