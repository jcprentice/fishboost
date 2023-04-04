# Full Pipeline ----
#
# Send simulated data to SIRE 2.x
# Eventually intend to replace simulated data with Turbot data
#
# SIRE 2.x modified from Chris Pooley's original

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
    row_no <- 1L
    protocol_file <- "param_sets/sim-patch-fb-1.rds"
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
            # "prior__abc_xyz__val1" -> "abc_xyz, val1"
            p1 <- unlist(strsplit(param, "__"))[-1]
            params$priors[parameter == p1[1], (p1[2]) := tmp]
        } else if (param == "group_effect") {
            params[param] <- tmp
            params$priors[parameter == "sigma", `:=`(
                val2 = max(2 * tmp, 1e-3),
                true_val = max(tmp, 0))]
        } else if (startsWith(param, "fe_")) {
            # fe_trial_infectivity -> "trial", "infectivity"
            p1 <- unlist(strsplit(param, "_"))[-1]
            td <- p1[1]; tr <- p1[2]; tr1 <- substr(p1[2], 1, 1)
            params$fe_vals[td, tr] <- tmp
            params$priors[parameter == paste0(td, "_", tr1), `:=`(
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
    params$nsample  <- as.integer(params$nsample)
    params$burnin   <- as.integer(params$burnin)
    params$thin     <- as.integer(params$thin)
    # overrule whatever is in the protocol and just take it from qsub
    params$nthreads <- as.integer(Sys.getenv("NUM_SLOTS"))
    
    # This can go awry
    if (is.na(params$group_effect)) {
        params$group_effect <- -1
    }
    
    # Patch params with posterior mean values from data set/scenario
    params <- patch_params(params)
    
    # fix any traits that need to be linked
    params <- apply_links(params)
    
    # Quick check of how things are looking
    summarise_params(params)
}


## Generate pedigree and traits ----

if (params$use_fb_data) {
    # Load pop, pedigree, and GRM
    pop <- switch(params$setup,
                  "fishboost" = readRDS("fb_data/fb_data12.rds"),
                  "fb1"       = readRDS("fb_data/fb_data1.rds"),
                  "fb2"       = readRDS("fb_data/fb_data2.rds"),
                  stop("Not a fishboost data set"))
    
    pedigree <- NULL
    GRM <- NULL
    
    # FB only has the last 2 events
    params$pass_events <- 2
} else {
    pedigree <- make_pedigree(params)
    
    # this is an A-matrix, not actually a GRM
    GRM <- NULL #make_grm(pedigree)
    
    # traits <- make_traits_from_grm(GRM, pedigree, params)
    if (params$patch_traits) {
        traits <- patch_in_traits(pedigree, params)
    } else {
        traits <- make_traits_from_pedigree(pedigree, params)
    }
    
    # set groups and donors
    traits <- set_groups(traits, params)
    
    # Set donor and trial effects
    traits <- apply_fixed_effects(traits, params)
}


## Run epidemic ----

if (params$use_fb_data) {
    # If we were plotting we might add missing timings here
    params$tmax <- get_tmax(pop, params)
} else {
    pop <- simulate_epidemic(traits, params)
    pop[sdp == "progeny", parasites := !is.na(Tinf)]
    params$estimated_R0 <- get_R0(pop)
    params$tmax <- get_tmax(pop, params)
}


## Generate XML File and run SIRE ----

{
    if (!dir.exists(params$data_dir)) {
        dir.create(params$data_dir)
    }
    
    # Tidy up ready for saving data to XML file
    data <- prepare_data(pop, params)
    # generate_sire2x_xml(data, params, GRM)
    get(glue("generate_{params$sire_version}_xml"))(data, params, GRM)
    
    cmd <- with(params, glue("../{sire_version}/sire {data_dir}/{name}.xml 0"))
    if (params$algorithm == "pas") {
        cmd <- glue("mpirun -n {params$nchains} --oversubscribe {cmd}")
    }
    
    # Run SIRE
    message(glue("Running {params$sire_version} ..."))
    time_taken <- system.time(system(cmd))
}


## Retrieve results and save ----

{
    message("Retrieving results ...")
    
    file_prefix <- with(params, glue("{data_dir}/{name}_out"))
    
    if (params$sire_version == "sire22") {
        parameter_estimates <- fread(glue("{file_prefix}/posterior.csv"))
        file_name           <- glue("{file_prefix}/ebvs.csv")
        estimated_BVs       <- if (file.exists(file_name)) fread(file_name) else NULL
        pred_accs           <- fread(glue("{file_prefix}/pred_accs.csv"))
    } else {
        parameter_estimates <- fread(glue("{file_prefix}/posteriors.csv"))
        estimated_BVs       <- fread(glue("{file_prefix}/estimated_bvs.csv"))
        pred_accs           <- fread(glue("{file_prefix}/pred_accs.csv"))
    }
    
    # ranks <- get_ranks(pop, estimated_BVs, params)
    message("parameter estimates:")
    msg_pars(parameter_estimates)
    
    if (!dir.exists(params$results_dir)) {
        dir.create(params$results_dir)
    }
    
    save(params, pop, parameter_estimates, estimated_BVs, pred_accs, time_taken, # ranks
         file = with(params, glue("{results_dir}/{name}.RData")))
}

message("Finished!")
