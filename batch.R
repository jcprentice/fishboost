# Full Pipeline ----
#
# Send simulated data to SIRE 2.0 / 2.1
# Eventually intend to replace simulated data with Turbot data
#
# SIRE 2.0 / 2.1 modified from Chris Pooley's original

## Load libraries and source files ----
source("libraries.R")
source("source_files.R")


# set.seed(0)

## Set parameters ----

cmd_args <- commandArgs(trailingOnly = TRUE)
if (length(cmd_args) == 2) {
    protocol_file <- glue("param_sets/{cmd_args[1]}.rds")
    row_no <- as.integer(cmd_args[2])
} else {
    row_no <- 3L
    protocol_file <- "param_sets/fb-mpi.rds"
}

{
    # grab row_no from protocol file
    message(glue("Running {protocol_file} / scenario {row_no}"))
    protocol <- readRDS(protocol_file)$protocol[row_no]

    # create params as usual
    params <- with(protocol,
                   make_parameters(
                       model_type = model_type,
                       data_set = data_set,
                       name = name,
                       use_traits = use_traits,
                       vars = vars,
                       covars = covars,
                       setup = setup,
                       group_effect = group_effect,
                       trial_fe = trial_fe,
                       donor_fe = donor_fe,
                       txd_fe = txd_fe))

    message("Copying parameters ...")

    # patch params from protocol file
    for (param in names(protocol)) {
        tmp <- protocol[[param]]

        # skip NAs
        if (is.na(tmp)) next

        message(" - ", param, " = ", tmp)

        # handle special cases first
        if (param == "Sigma") {
            tmp <- eval(parse(text = tmp))
            params$Sigma_G[,] <- tmp
            params$Sigma_E <- params$Sigma_G
        } else if (param == "latency_period") {
            params$latency_period <- tmp
            params$r_eta <- 1.0 / tmp
            params$r_eta_shape <- 1.0 # 10.0
            params$r_eta_rate <- params$r_eta * params$r_eta_shape
        } else if (param == "detection_period") {
            params$detection_period <- tmp
            params$r_rho <- 1.0 / tmp
            params$r_rho_shape <- 1.0 # 10.0
            params$r_rho_rate <- params$r_rho * params$r_rho_shape
        } else if (param == "recovery_period") {
            params$recovery_period <- tmp
            params$r_gamma <- 1.0 / tmp
            params$r_gamma_shape <- 1.0 # 10.0
            params$r_gamma_rate <- params$r_gamma * params$r_gamma_shape
        } else if (startsWith(param, "prior")) {
            # "prior__abc_xyz__1" -> "abc_xyz"
            p1 <- unlist(strsplit(param, "__"))[2:3]
            params$priors[parameter == p1[1], paste0("val", p1[2]) := tmp]
        } else if (param == "group_effect") {
            params[param] <- tmp
            params$priors[parameter == "sigma", `:=`(
                val2 = max(2 * tmp, 1e-3),
                true_val = max(tmp, 0))]
        } else if (startsWith(param, "fe_")) {
            # fe_trial_infectivity -> "trial", "infectivity"
            p1 <- unlist(strsplit(param, "_"))[2:3]
            td <- p1[1]; tr <- p1[2]
            params$fe_vals[td, tr] <- tmp
            params$priors[parameter == paste0(td, "_", substr(tr, 1, 1)), `:=`(
                val1 = min(tmp - 2, val1),
                val2 = max(tmp + 2, val2),
                true_val = tmp)]
        } else if (endsWith(param, "_fe")) {
            # don't accidentally reset this to "none"
            params[param] <- if (tmp == "none") "" else tmp
        } else if (param == "seed" && tmp > 0) {
            set.seed(tmp)
        } else {
            params[param] <- tmp
        }
    }

    # tidy up numbers (only in balanced / non-Fishboost cases)
    if (!params$setup %in% c("fishboost", "fb1", "fb2")) {
        params$ndams      <- with(params, nsires * dpsire)
        params$nparents   <- with(params, nsires + ndams)
        params$nprogeny   <- with(params, ndams * ppdam)
        params$ntotal     <- with(params, nparents + nprogeny)
        params$group_size <- with(params, nprogeny / ngroups)
    }
    params$R0 <- with(params, r_beta / r_gamma)

    # ensure these are integers!
    params$nsample <- as.integer(params$nsample)
    params$burnin  <- as.integer(params$burnin)
    params$thin    <- as.integer(params$thin)
    
    # This can go awry
    if (is.na(params$group_effect)) {
        params$group_effect <- -1
    }

    # fix any traits that need to be linked
    params <- apply_links(params)
}

summarise_params(params)

## Generate pedigree and traits ----

if (params$use_fb_data) {
    # Load pop, pedigree, and GRM
    switch(params$setup,
           "fishboost" = load("fb_data/fb_traits.RData"),
           "fb1" = load("fb_data/fb_traits1.RData"),
           "fb2" = load("fb_data/fb_traits2.RData"),
           stop("Not a fishboost data set"))
    
    # FB only has the last 2 events
    params$pass_events <- 2
} else {
    pedigree <- make_pedigree(params)

    # this is an A-matrix, not actually a GRM
    GRM <- make_grm(pedigree)

    # traits <- make_traits_from_grm(GRM, pedigree, params)
    traits <- make_traits_from_pedigree(pedigree, params)

    # set groups and donors
    traits <- set_groups(traits, params)

    # Set donor and trial effects
    traits <- apply_fixed_effects(traits, params)
}


## Run epidemic ----

if (params$use_fb_data) {
    # If we were plotting we might add missing timings here
} else {
    pop <- simulate_epidemic(traits, params)
    params$estimated_R0 <- get_R0(pop)
}


## Generate XML File and run SIRE ----

{
    if (!dir.exists(params$data_dir)) {
        dir.create(params$data_dir)
    }

    # Tidy up ready for saving data to XML file
    data <- prepare_data(pop, params)

    # Generate the XML file
    if (params$sire_version == "2.2") {
        x_xml <- generate_sire22_xml(data, GRM, params)
        # do we need --oversubscribe here?1
        cmd <- with(params, glue("mpirun -n {nthreads} --oversubscribe ../sire22/sire {data_dir}/{name}.xml {replicate}"))
        message("Running SIRE 2.2 ...")
    } else {
        x_xml <- generate_sire21_xml(data, GRM, params)
        cmd <- with(params, glue("../sire21/sire {data_dir}/{name}.xml {replicate}"))
        message("Running SIRE 2.1 ...")
    }

    # Run SIRE
    time_taken <- system.time(system(cmd))
}


## Retrieve results and save ----

{
    message("Retrieving results ...")

    # note: sire 2.0 needs "data_dir/name-", instead of "name_out/"
    file_prefix <- with(params, glue("{data_dir}/{name}_out"))

    if (params$sire_version == "2.2") {
        parameter_estimates <- fread(glue("{file_prefix}/posterior.csv"))
        file_name <- glue("{file_prefix}/ebvs.csv")
        estimated_BVs       <- if (file.exists(file_name)) fread(file_name) else NULL
        pred_accs           <- fread(glue("{file_prefix}/pred_accs.csv"))
    } else {
        parameter_estimates <- fread(glue("{file_prefix}/posteriors.csv"))
        estimated_BVs       <- fread(glue("{file_prefix}/estimated_bvs.csv"))
        pred_accs           <- fread(glue("{file_prefix}/pred_accs.csv"))
    }

    # ranks <- get_ranks(pop, estimated_BVs, params)

    if (!dir.exists(params$results_dir)) {
        dir.create(params$results_dir)
    }

    save(params, pop, parameter_estimates, estimated_BVs, pred_accs, time_taken, # ranks
         file = with(params, glue("{results_dir}/{name}.RData")))
}

message("Finished!")
