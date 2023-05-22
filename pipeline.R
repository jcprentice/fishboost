# Full Pipeline ----
#
# Either use the Fishboost data set or simulate a new data set, and then send it
# to SIRE 2.2 for parameter inference.
#
# SIRE 2.2 (modified from Chris Pooley's original)

## Load libraries and source files ----

source("libraries.R")
source("source_files.R")

# for testing only
# set.seed(0)


## Set parameters ----

{
    params <- make_parameters(
        model_type = "SEIDR", # "SIR", "SEIR", "SIDR", or "SEIDR"
        name = "expt1",
        setup = "fishboost", # chris, small, fishboost, fb1, fb2
        use_traits = "sir", # "all", "none", "sir", "si" etc.
        vars = 0.5, # c(0.5, 1.2, 0.5)
        covars = 0.2, # 0.2, "positive", "negative", "mixed", "sir_only" etc.
        group_layout = "fishboost", # "random", "family", "striped", "fishboost"
        trial_fe = "lidr",
        donor_fe = "lidr",
        # txd_fe = "lidr",
        use_fb_data = TRUE
    )
    
    # Temporary override of some parameters
    # params$group_effect <- 0.3
    params$nthreads <- 2L
    params$nsample <- 1e4L
    params$burnin <- 2e3L
    params$thin <- 1e0L
    params$nsample_per_gen <- 1e1L
    params$anneal <- "on"
    params$anneal_power <- 4
    params$pass_events <- 2
    # params$link_traits <- "srirr"
    # params$sim_link_donor <- "slill"
    # params$link_donor <- "slill"
    params$time_step <- 1
    params$censor <- 0.8
    
    # Patch params with posterior mean values from data set/scenario
    params$patch_data_set <- "fb-parasites"
    params$patch_scenario <- 1
    params <- patch_params(params)

    # link traits (needs to be done after params is constructed)
    params <- apply_links(params)
    
    # Quick check of how things are looking
    summarise_params(params)
}


## Generate pedigree and traits ----

if (params$use_fb_data) {
    fb_traits <- switch(params$setup,
                        "fishboost" = readRDS("fb_data/fb_data12.rds"),
                        "fb1"       = readRDS("fb_data/fb_data1.rds"),
                        "fb2"       = readRDS("fb_data/fb_data2.rds"))
    
    # fb_traits[Tsym == 1, Tsym := 2]
    
    # Create a GRM
    GRM <- NULL #make_grm(fb_traits)
    
    # We don't know the BVs, so just set them to whatever (SIRE doesn't seem to like NA)
    # fb_traits[, `:=`(susceptibility_BV = NA_real_,
    #                  infectivity_BV    = NA_real_,
    #                  recoverability_BV = NA_real_)]
    
} else {
    pedigree <- make_pedigree(params)
    
    # This is an A-matrix, not actually a GRM
    GRM <- NULL #make_grm(pedigree)
    
    # traits <- make_traits_from_grm(GRM, pedigree, params)
    if (params$patch_traits) {
        traits <- patch_in_traits(pedigree, params)
    } else {
        traits <- make_traits_from_pedigree(pedigree, params)
    }
    
    # Set groups, trial, donors, and group effect
    traits <- set_groups(traits, params)
    
    # Set donor and trial effects
    traits <- apply_fixed_effects(traits, params)
}


## Simulate and plot epidemic ----

if (params$use_fb_data) {
    pop <- fb_traits
    
    # If Tinf, Tinc etc. are not in pop, create missing column(s) and set to NA
    missing_timings <- setdiff(params$timings, names(pop))
    if (length(missing_timings) > 0) {
        pop[, (missing_timings) := NA_real_]
    }
    
    params$tmax <- get_tmax(pop, params)
    pop[Trec == params$tmax[trial], Trec := NA]
    
    plt <- plot_model(pop, params)
    # plt_fb1 <- plot_SxxDR(pop, params)
} else {
    pop <- simulate_epidemic(traits, params)
    pop[sdp == "progeny", parasites := !is.na(Tinf)]
    
    params$estimated_R0 <- get_R0(pop)
    params$tmax <- get_tmax(pop, params)
    
    ### Plot the epidemic ----
    plt <- plot_model(pop, params)
}


# Generate XML File and run SIRE ----

{
    if (!dir.exists(params$data_dir)) {
        dir.create(params$data_dir)
    }
    
    # Tidy up ready for saving data to XML file
    params$use_parasites <- "SE" # S, SE, or SEI
    params$time_step
    data3 <- prepare_data(pop, params)
    
    # generate_sire2x_xml(data, params, GRM)
    sire_xml <- get(glue("generate_{params$sire_version}_xml"))(data, params, GRM)
    
    cmd <- with(params, glue("../{sire_version}/sire {data_dir}/{name}.xml 0"))
    if (params$algorithm == "pas") {
        cmd <- glue("mpirun -n {params$nchains} --oversubscribe {cmd}")
    }
    
    # Run SIRE
    message(glue("Running: {cmd} ..."))
    time_taken <- system.time(system(cmd))
}


## Retrieve results and save ----

{
    message("Retrieving results ...")
    
    # note: sire 2.0 needs "dat_dir/name-", instead of "name_out/"
    file_prefix <- with(params, glue("{data_dir}/{name}_out"))
    
    if (params$sire_version == "2.2") {
        parameter_estimates <- fread(glue("{file_prefix}/posterior.csv")) # 2.1 uses posteriors.csv
        file_name <- glue("{file_prefix}/ebvs.csv") # 2.1 uses estimated_bvs.csv
        estimated_BVs       <- if (file.exists(file_name)) fread(file_name) else NULL
        pred_accs           <- fread(glue("{file_prefix}/pred_accs.csv"))
    } else {
        parameter_estimates <- fread(glue("{file_prefix}/posteriors.csv"))
        estimated_BVs       <- fread(glue("{file_prefix}/estimated_bvs.csv"))
        pred_accs           <- fread(glue("{file_prefix}/pred_accs.csv"))
    }
    
    msg_pars(parameter_estimates)
    print(parameter_estimates[!startsWith(parameter, "Group effect")])
    if (!params$use_fb_data) {
        print(pred_accs)
    }
    
    ranks <- get_ranks(pop, estimated_BVs, params)
    
    ### Save results ----
    if (!dir.exists(params$results_dir)) {
        dir.create(params$results_dir)
    }
    
    save(params, pop, parameter_estimates, estimated_BVs, pred_accs, ranks, time_taken,
         file = with(params, glue("{results_dir}/{name}.RData")))
}
