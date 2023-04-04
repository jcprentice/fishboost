library(data.table)

# Generate parameter list
make_parameters <- function(
    # Epidemic model: "SIR", "SIDR", "SEIR", "SEIDR", "...-res" for reservoir model
    model_type = "SEIDR",
    # Filename for data
    name = "expt",
    # Folder suffix for data (default "")
    data_set = "",
    # Population layout
    setup = "fishboost",
    # Which traits to use (this is going to be "clever")
    # use "all" or "none", or a subset using the first letter of each trait
    # "slidr" = "sus", "lat", "inf", "det", "rec"
    use_traits = "sir",
    # Variances: diagonal of cov matrix, either a single value, or a vector, in
    # which case it should be compatible with "use_traits"
    vars = 1,
    # Cov matrix, numerical: 0.2, c(0.2, 0.2), or string: "chris", "positive", "sir_only"
    covars = 0.2,
    # How to arrange individuals into groups: "random", "striped", "family", "fishboost"
    group_layout = "fishboost",
    # Include group effect: if >= 0 simulate if appropriate and check, ignore if < 0
    group_effect = -1,
    # List of fixed effects for SIRE to check for donors and trial. Should also
    # simulate these if not using FB data_set A string with the first letters of traits.
    trial_fe = "none",
    donor_fe = "none",
    txd_fe = "none",
    # Simulate new data or use Fishboost data
    use_fb_data = FALSE
) {

    message("Setting parameters ...")

    # Fix input NAs ----

    # Handy in case of debugging
    if (FALSE) {
        model_type <- "SEIDR"; name <- "expt"; data_set <- ""; setup <- "fishboost";
        use_traits <- "si"; vars <- 1; covars <- 0.2; group_layout <- "random";
        group_effect <- -1; trial_fe <- "lid"; donor_fe <- "lid"; txd_fe <- "lid";
        use_fb_data <- FALSE;
    }
    
    # Handy for batch
    if (FALSE) {
        model_type <- protocol$model; data_set <- protocol$data_set; name <- protocol$name;
        use_traits <- protocol$use_traits; vars <- protocol$vars; covars <- protocol$covars;
        setup <- protocol$setup; group_layout <- protocol$group_layout; group_effect <- protocol$group_effect;
        donor_fe <- protocol$donor_fe; trial_fe <- protocol$trial_fe; txd_fe <- protocol$txd_fe;
        use_fb_data <- protocol$use_fb_data
    }

    # Check if value is unassigned, NA, or NULL, and replace with sensible defaults
    get_default <- function(x, x_def) {
        if (is.na(x) || is.null(x)) x_def else x
    }

    model_type   <- get_default(model_type, "SIR")
    name         <- get_default(name, paste0(tolower(model_type), 1L))
    data_set     <- get_default(data_set, "")
    setup        <- get_default(setup, "none")
    use_traits   <- get_default(use_traits, "sir")
    vars         <- get_default(vars, 1)
    covars       <- get_default(covars, 0.5)
    group_layout <- get_default(group_layout, "fishboost")
    trial_fe     <- get_default(trial_fe, "none")
    donor_fe     <- get_default(donor_fe, "none")
    txd_fe       <- get_default(txd_fe, "none")
    use_fb_data  <- get_default(use_fb_data, FALSE)

    # Set output directories ----

    # Directories should be something like "data", "data-sim-xyz/", or
    # "data-fb-xyz/", so make sure that "-" is in there
    data_dir    <- paste0("data", if (data_set != "") "-" else "", data_set)
    results_dir <- sub("data", "results", data_dir)


    # Population setup ----
    message(" - Setup is: '", setup, "'")
    message(" - Group layout is: '", group_layout, "'")

    # Note: in balanced populations, dpsire = dams per sire, ppdam = progeny per dam
    switch(setup,
           "fishboost" = {
               nsires <- 29L; ndams <- 25L; nprogeny <- 1775L;
               dpsire <- 1L; ppdam <- 72L;
               ngroups <- 71L; I0 <- 5L;
           }, "fb1" = {
               nsires <- 14L; ndams <- 14L; nprogeny <- 875L;
               dpsire <- 1L; ppdam <- 72L;
               ngroups <- 35L; I0 <- 5L;
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

    # Derived numbers
    nparents <- nsires + ndams
    ntotal <- nprogeny + nparents
    group_size <- nprogeny / ngroups


    # Traits ----
    message(" - Model type is: '", model_type, "'")

    # traitnames        list of traits that the model uses
    # corr_signs        for when covariance is mixed, should more disease = higher (+1) or lower (-1)
    # compartments      names of compartments
    # timings           name for time when individual enters compartment
    # use_parameters    subset of parameters to pass to SIRE 2.1

    # Check for reservoir model
    reservoir <- endsWith(tolower(model_type), "-res")
    if (reservoir) {
        model_type <- sub("-res", "", model_type)
    }

    # Compartments are just the letters in model_type (ignore repeated S)
    compartments <- unique(strsplit(model_type, "")[[1]])

    complete_traitnames <- c("susceptibility", "latency", "infectivity", "detectability", "recoverability")

    # Set up which traits a model will use, if traits are +vely or -vely
    # correlated, what the event times are called, and which parameters we need
    # to make SIRE aware of.
    switch(
        model_type,
        "SEIDR" = {
            all_traitnames <- c("susceptibility", "latency", "infectivity", "detectability", "recoverability")
            corr_signs <- c(1, 1, -1, 1, 1)
            # corr_signs <- c(1, -1, 1, -1, -1)
            timings <- c("Tinf", "Tinc", "Tsym", "Trec")
            use_parameters <- c("beta", "latent_period", "eta_shape", "detection_period",
                                "rho_shape", "recovery_period", "gamma_shape")
        }, "SIDR" = {
            all_traitnames <- s("susceptibility", "infectivity", "detectability", "recoverability")
            corr_signs <- c(1, 1, -1, -1)
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
    
    # Sigma_G <- matrix(0, 5, 5, dimnames = list(all_traitnames, all_traitnames))

    # Likely using only a subset of traits
    traitnames <- if (use_traits == "all") {
        all_traitnames
    } else if (use_traits == "none") {
        NULL
    } else {
        all_traitnames[match(strsplit(use_traits, "")[[1]], substr(all_traitnames, 1, 1))]
    }
    
    
    names(corr_signs) <- all_traitnames
    corr_signs <- unname(corr_signs[traitnames])
    
    ntraits <- length(traitnames)
    nall_traits <- length(all_traitnames)

    message(" - ", ntraits, " traits: ", paste0(traitnames, collapse = ", "))
    
    # Assume default of var(i) = 0.5 and cov(i, j) = 0
    Sigma_G <- diag(vars[vars > 0], nall_traits, nall_traits)
    dimnames(Sigma_G) <- list(all_traitnames, all_traitnames)
    Sigma_E <- Sigma_G


    # Useful to have these indices (since they depend on the model). I could
    # just write out the trait name each time, but this is easier!
    iS <- which(all_traitnames == "susceptibility")
    iI <- which(all_traitnames == "infectivity")
    iR <- which(all_traitnames == "recoverability")

    iSIR <- c(iS, iI, iR)
    iSI  <- c(iS, iI)
    
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
                Sigma_G[iSIR, iSIR] <- c(+0.33,  +0.3,  -0.044,
                                         +0.3,   +1.68, -0.5,
                                         -0.044, -0.5,  +0.6);
                Sigma_E[iSIR, iSIR] <- c(+0.77,  +0.56, -0.25,
                                         +0.56,  +1.12, -0.4,
                                         -0.25,  -0.4,  +0.9)
            }, "positive" = {
                Sigma_G[] <- 0.4
                diag(Sigma_G) <- 1.0
            }, "negative" = {
                Sigma_G[] <- -0.4
                diag(Sigma_G) <- 1.0
            }, "mixed" = {
                Sigma_G[] <- 0.4
                diag(Sigma_G) <- 1.0
                Sigma_G <- Sigma_G * corr_signs %*% t(corr_signs)
            }, "sir_weak" = {
                Sigma_G[iSIR, iSIR] <- 0.1
                diag(Sigma_G) <- 1e-6
                diag(Sigma_G[iSIR, iSIR]) <- 1.0
                Sigma_G <- Sigma_G * corr_signs %*% t(corr_signs)
            }, "sir_no_cor" = {
                diag(Sigma_G) <- 1e-6
                diag(Sigma_G[iSIR, iSIR]) <- 1.0
            }, "sir_pos_cor" = {
                Sigma_G[iSIR, iSIR] <- 0.4
                diag(Sigma_G) <- 1e-6
                diag(Sigma_G[iSIR, iSIR]) <- 1.0
            }, "si_no_cor" = {
                diag(Sigma_G) <- 1e-6
                diag(Sigma_G[iSI, iSI]) <- 1.0
            }, "si_pos_cor" = {
                Sigma_G[iSI, iSI] <- 0.4
                diag(Sigma_G) <- 1e-6
                diag(Sigma_G[iSI, iSI]) <- 1.0
            }, "sir_only" = {
                diag(Sigma_G) <- 1e-6
                diag(Sigma_G[iSIR, iSIR]) <- 1.0
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
    cor_G <- cov2cor(Sigma_G)
    cor_E <- cov2cor(Sigma_E)


    # Epidemic parameters ----

    ## Infection coefficient ----
    r_beta <- 0.2

    ## Reservoir ----
    r_zeta   <- 0.05 # shedding
    r_lambda <- 0.5 # decay

    # Note for Gamma dist:
    # shape = mean^2 / var
    # rate  = mean / var

    ## Latent Periods ----

        latent_period <-  2

    ### Incubation period ----

    # Calculated for an SEIDR model based on donors in trial 1 (+ 0.5)
    # fitdist(Tsym + 0.5, "gamma")
    r_eta <- 1 / latent_period
    r_eta_shape <- 1 # 1.46 # 2.0 # 10.0
    r_eta_rate  <- r_eta * r_eta_shape # 0.232 # r_eta_shape / latent_period
    eta_type <- "exp"


    ### Asymptomatic infectious period ----

    # This is currently identical to the latency
    r_rho <- 1 / latent_period
    r_rho_shape <- 1 # 1.46 # 2.0 # 10.0
    r_rho_rate  <- r_rho * r_rho_shape # 0.232 # r_eta_shape / latent_period
    rho_type <- "exp"


    ## Recovery ----

    # fitdist(RP + 0.5, "gamma")
    recovery_period <- 8
    r_gamma <- 1 / recovery_period
    r_gamma_shape <- 1
    r_gamma_rate <- r_gamma_shape / recovery_period # 0.149
    gamma_type <- "exp"


    ## Fixed effects ----

    # prevent fread interpreting these as NA
    if (trial_fe == "none") trial_fe <- ""
    if (donor_fe == "none") donor_fe <- ""
    if (txd_fe   == "none") txd_fe   <- ""
    
    # when simulating, copy these values, or override them after params is
    # created if they should be different
    sim_trial_fe <- trial_fe
    sim_donor_fe <- donor_fe
    sim_txd_fe   <- txd_fe
    
    # tell SIRE to include these variables
    use_parameters <- c(
        use_parameters,
        if (trial_fe == "") NULL else paste0("trial_", strsplit(trial_fe, "")[[1]]),
        if (donor_fe == "") NULL else paste0("donor_", strsplit(donor_fe, "")[[1]]),
        if (txd_fe   == "") NULL else paste0("txd_",   strsplit(txd_fe,   "")[[1]])
    )

    # c(0, 3,     1.5, 1, 0,
    #   0, 1.75, -1,   1, 0)
    #                   s   l   i  d  r
    fe_vals <- matrix(c(0,  1, -1, 1, 0,  # trial
                        0, -1,  1, 1, 0,  # donor
                        0,  1,  1, 1, 0), # txd
                        # 0, -3, 2, 3, 2), # donor
                      nrow = 3, ncol = 5, byrow = TRUE,
                      dimnames = list(c("trial", "donor", "txd"),
                                      complete_traitnames))

    ## Map traits and FEs ----
    
    # "xxxxx" will be copied onto "slidr", so "srirr" means that lat, det, rec
    # will all use rec. sim means simulate with, otherwise it's what SIRE sees
    sim_link_traits <- "slidr"
    sim_link_trial <- "slidr"
    sim_link_donor <- "slidr"
    sim_link_txd <- "slidr"

    link_traits <- "slidr"
    link_trial <- "slidr"
    link_donor <- "slidr"
    link_txd <- "slidr"
    
    # Shape parameters
    sim_link_shapes <- "ldr"
    link_shapes <- "ldr"
    
    ## R0 ----

    # This should target R0 around 5?
    R0 <- r_beta / (r_gamma + if ("D" %in% compartments) r_rho else 0)

    # MCMC settings ----
    
    # Which version of SIRE TO USE
    sire_version <- "sire22"
    
    nchains <- 10L

    ## Samples ----
    nsample <- as.integer(1e6)
    burnin  <- as.integer(nsample / 5)
    thin    <- as.integer(max(nsample / 1e4, 1))

    nsample_per_gen <- 1e3L

    nthreads <- 1L
    
    use_pas <- TRUE
    if (use_pas) {
        algorithm <- "pas"; anneal <- NULL; anneal_power <- NULL
    } else {
        algorithm <- "mcmc"; anneal <- "on"; anneal_power <- 4
    }
    
    phi <- 1.0

    ## Priors ----
    
    if (FALSE) {
        keep_traits <- all_traitnames[all_traitnames %in% traitnames]
        drop_traits <- all_traitnames[!all_traitnames %in% traitnames]
        
        Sigma_G[drop_traits, drop_traits] <- NA
        Sigma_E[drop_traits, drop_traits] <- NA
        cor_G[drop_traits, drop_traits] <- NA
        cor_E[drop_traits, drop_traits] <- NA
    }
    
    priors <- as.data.table(dplyr::tribble(
        ~parameter,         ~type,  ~val1, ~val2, ~true_val,       ~use,
        
        "beta",             "Flat", 0,     1.0,   r_beta,          TRUE,
        "latent_period",    "Flat", 0,     160,   latent_period,   FALSE,
        "eta_shape",        "Flat", 0.99,  1.01,  r_eta_shape,     FALSE,
        "detection_period", "Flat", 0,     160,   latent_period,   FALSE,
        "rho_shape",        "Flat", 0.99,  1.01,  r_rho_shape,     FALSE,
        "recovery_period",  "Flat", 0,     20,    recovery_period, FALSE,
        "gamma_shape",      "Flat", 0.5,   5,     r_gamma_shape,   FALSE,
        "sigma",            "Flat", 0,     0.3,   group_effect,    TRUE,
        
        # parameter  type    val1  val2  true_val         use
        "cov_G_ss", "Flat",  1e-3, 5,    Sigma_G[iS, iS], TRUE,
        "cov_G_ii", "Flat",  1e-3, 4,    Sigma_G[iI, iI], TRUE,
        "cov_G_rr", "Flat",  1e-3, 1,    Sigma_G[iR, iR], TRUE,
        "r_G_si",   "Flat", -0.95, 0.95, cor_G[iS, iI],   TRUE,
        "r_G_sr",   "Flat", -0.95, 0.95, cor_G[iS, iR],   TRUE,
        "r_G_ir",   "Flat", -0.95, 0.95, cor_G[iI, iR],   TRUE,
        "cov_E_ss", "Flat",  1e-3, 5,    Sigma_E[iS, iS], TRUE,
        "cov_E_ii", "Flat",  1e-3, 4,    Sigma_E[iI, iI], TRUE,
        "cov_E_rr", "Flat",  1e-3, 1,    Sigma_E[iR, iR], TRUE,
        "r_E_si",   "Flat", -0.95, 0.95, cor_E[iS, iI],   TRUE,
        "r_E_sr",   "Flat", -0.95, 0.95, cor_E[iS, iR],   TRUE,
        "r_E_ir",   "Flat", -0.95, 0.95, cor_E[iI, iR],   TRUE,
        
        # parameter type   val1 val2 true_val                            use
        "trial_l", "Flat", -4,  8,   fe_vals["trial", "latency"],        grepl("l", trial_fe),
        "trial_i", "Flat", -3,  5,   fe_vals["trial", "infectivity"],    grepl("i", trial_fe),
        "trial_d", "Flat", -2,  7,   fe_vals["trial", "detectability"],  grepl("d", trial_fe),
        "trial_r", "Flat", -1,  1,   fe_vals["trial", "recoverability"], grepl("r", trial_fe),
        "donor_l", "Flat", -6,  6,   fe_vals["donor", "latency"],        grepl("l", donor_fe),
        "donor_i", "Flat", -4,  4,   fe_vals["donor", "infectivity"],    grepl("i", donor_fe),
        "donor_d", "Flat", -4,  8,   fe_vals["donor", "detectability"],  grepl("d", donor_fe),
        "donor_r", "Flat", -4,  2,   fe_vals["donor", "recoverability"], grepl("r", donor_fe),
        "txd_l",   "Flat", -4, 10,   fe_vals["txd",   "latency"],        grepl("l", txd_fe),
        "txd_i",   "Flat", -6,  6,   fe_vals["txd",   "infectivity"],    grepl("i", txd_fe),
        "txd_d",   "Flat", -4,  12,  fe_vals["txd",   "detectability"],  grepl("d", txd_fe),
        "txd_r",   "Flat",  0,  4,   fe_vals["txd",   "recoverability"], grepl("r", txd_fe)
    ))
    
    # Drop unwanted traits from cov matrices
    if (FALSE) {
        Sigma_G <- Sigma_G[keep_traits, keep_traits, drop = FALSE]
        Sigma_E <- Sigma_E[keep_traits, keep_traits, drop = FALSE]
        cor_G <- cor_G[keep_traits, keep_traits, drop = FALSE]
        cor_E <- cor_E[keep_traits, keep_traits, drop = FALSE]
        h2 <- h2[keep_traits, keep_traits, drop = FALSE]
    }
    
    # Some priors set to not use with SIRE, update if we really do want them
    
    # # First remove exp / gamma shaped periods
    # cut_types <- c(if (eta_type == "exp") "eta_shape" else NULL,
    #                if (rho_type == "exp") "rho_shape" else NULL,
    #                if (gamma_type == "exp") "gamma_shape" else NULL)
    # use_parameters <- use_parameters[!use_parameters %in% cut_types]

    # Then set all use values
    priors[parameter %in% use_parameters, use := TRUE]
    

    # remove priors if traits not used
    if (!"susceptibility" %in% traitnames) {
        priors[parameter %in% c("cov_G_ss", "r_G_si", "r_G_sr", "cov_E_ss", "r_E_si", "r_E_sr"), use := FALSE]
    }
    if (!"infectivity" %in% traitnames) {
        priors[parameter %in% c("cov_G_ii", "r_G_si", "r_G_ir", "cov_E_ii", "r_E_si", "r_E_ir"), use := FALSE]
    }
    if (!"recoverability" %in% traitnames) {
        priors[parameter %in% c("cov_G_rr", "r_G_sr", "r_G_ir", "cov_E_rr", "r_E_sr", "r_E_ir"), use := FALSE]
    }

    # Additional settings ----

    # What we pass to SIRE 2.0, choice of
    # "Tinf": actual infection time
    # "Tsym": time of symptoms
    # "estimated_Tinf_from_donors": Tsym - mean incubation period of donors
    # "estimated_Tinf_per_individual": Tsym with all gaps covered
    pass_Tsym <- "Tsym"
    
    # pass event times, pass the last N of timings (Tinf, Tinc, Tdet, Trec)
    pass_events <- 2 # pass time Tdet and Trec
    
    time_step <- 0

    censor <- 1.0
    
    # how to handle the parasites column, if parasites == F, then restrict
    # "S|E|I" to "S|E" or even "S". If parasites == T, but no symptoms, then
    # restrict to "S|E|I", "E|I", or "I".
    use_parasites <- ""

    # Patch data with values from from data_set / scenario posterior
    patch_data_set <- ""
    patch_scenario <- 1
    # Use means / randomly sample from posterior
    patch_with_mean <- FALSE
    # Patch traits with posterior EBVs
    patch_traits <- FALSE
    
    # Show plots during run
    show_plots <- TRUE

    # Show details
    DEBUG <- FALSE


    # Create params list ----
    #
    # This will be a list for passing to functions
    params <- mget(
        # inputs
        c("data_set", "name", "data_dir", "results_dir", "model_type",
          "use_fb_data", "setup", "vars", "covars", "group_layout", "use_traits",
          # population
          "nsires", "ndams", "nparents", "nprogeny", "ngroups", "ntotal",
          "dpsire", "ppdam", "group_size", "I0",
          # model traits
          "compartments", "reservoir", "timings", "all_traitnames", "traitnames",
          "ntraits", "Sigma_E", "Sigma_G", "cor_G", "cor_E", "h2",
          "sim_link_traits", "sim_link_trial", "sim_link_donor", "sim_link_txd", "sim_link_shapes",
          "link_traits", "link_trial", "link_donor", "link_txd", "link_shapes",
          # model parameters (rates)
          "r_beta", "r_zeta", "r_lambda",                                              # infection
          "latent_period", "r_eta", "r_eta_shape", "r_eta_rate", "eta_type",           # latency
          "r_rho", "r_rho_shape", "r_rho_rate", "rho_type",                            # detection
          "recovery_period", "r_gamma", "r_gamma_shape", "r_gamma_rate", "gamma_type", # recovery
          "R0", "group_effect", "fe_vals",
          "sim_trial_fe", "sim_donor_fe", "sim_txd_fe", "trial_fe", "donor_fe", "txd_fe",
          # mcmc & extra
          "sire_version", "algorithm", "anneal", "anneal_power", "nchains",
          "nsample", "burnin", "thin", "nthreads", "phi", "nsample_per_gen",
          "priors", "time_step", "censor", "use_parasites", "pass_events",
          "patch_data_set", "patch_scenario", "patch_with_mean", "patch_traits",
          "show_plots", "DEBUG"))

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
            "     R0 = {signif(R0, 5)}\n",
            " - Running MCMC with:\n",
            "     {nsample} samples / {burnin} burnin / {thin} thin / {nthreads} threads\n",
            " - Data will be saved to '{data_dir}/{name}.xml'"
        ))
    })
}
