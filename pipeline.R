# Full Pipeline ----
#
# Either use the Fishboost data set or simulate a new data set, and then send it
# to SIRE or BICI (Chris Pooley) for parameter inference.

## Load libraries and source files ----

suppressPackageStartupMessages({
    source("libraries.R")
    source("source_files.R")
})

time_start <- now()

# for testing only
# set.seed(0)

cmd_args <- commandArgs(trailingOnly = TRUE)
run_from_script <- length(cmd_args) > 0

## Set parameters ----

{
    params <- make_parameters(
        model_type = "SEIDR", # "SIR", "SEIR", "SIDR", or "SEIDR"
        dataset = "testing",
        name = "scen-1-1",
        setup = "fb_12_rpw", # chris, small, fb_12, fb_1, fb_2, single
        use_traits = "sit", # "all", "none", "sit", "si" etc.
        vars = 0.5, # c(0.5, 1.5, 0, 0, 0.5)
        cors = 0.2, # c(0.2, ..., 0.2)
        group_layout = "fishboost", # "random", "family", "striped", "fishboost"
        trial_fe = "ildt",
        donor_fe = "ildt",
        txd_fe = "ildt",
        weight_fe = "sildt",
        weight_is_nested = TRUE,
        sim_new_data = "r"
    )
    
    # Temporary override of some parameters
    # params$group_effect <- 0.3
    # params$priors[parameter == "latent_period", `:=`(type = "Fixed", val1 = 5, true_val = 5)]
    
    if (!run_from_script) {
        params$nchains <- 2L
        params$nsample <- 1e2L
    }
    params$sample_states <- 100L
    params$ie_output <- "true"
    params$pass_events <- c("Tsym", "Tdeath")
    # params$link_traits <- "sittt"
    # params$sim_link_donor <- "sittt"
    # params$link_donor <- "sittt"
    params$time_step <- 1
    params$censor <- 0.8
    params$use_weight <- "log"
    params$use_grm <- "pedigree"
    params$weight_is_nested <- TRUE
    # params$fix_donors <- c(params$fix_donors, "set_to_R")
    
    # Tidy up LP, DP, RP, including rate, shape, and priors
    params <- tidy_up_periods(params)
    
    # Patch params with posterior mean values from data set/scenario
    params$patch_dataset <- "" # "fb-test"
    params$patch_name <- "scen-1-1"
    params$traits_source <- "posterior"
    params$patch_type <- "median"
    params$patch_state <- TRUE
    params$skip_patches <- c("") #, "covariance")
    
    # params$ge_opts <- "e1" # Optional genetic effects
    
    params <- params |>
        patch_params() |>  # Patch params with posteriors from dataset / scenario
        # set_ge_opts() |>   # Genetic covariance options
        set_use_flags() |> # Ensure priors are correctly enabled
        apply_links()      # Fix any traits that need to be linked
    
    # Quick check of how things are looking
    summarise_params(params)
}


## Generate pedigree and popn ----

if (params$sim_new_data != "no") {
    popn <- make_pedigree(params) |>
        set_groups(params) |> # Set groups, trial, donors, and group effect
        set_traits(params) |>
        set_weights(params) |>
        apply_fixed_effects(params)
} else {
    popn <- readRDS(str_glue("fb_data/{params$setup}.rds"))
    
    params$fix_donors <- "no_Tsym_survivors" # c("time", "no_Tsym_survivors")
    
    # Create a GRM or A matrix
    GRM <- make_grm(popn, params$use_grm)
}


## Simulate and plot epidemic ----

if (params$sim_new_data == "r") {
    tic(); popn <- simulate_epidemic(popn, params); toc()
    popn[sdp == "progeny", parasites := !is.na(Tinf)]
    
    params$estimated_R0 <- get_R0(popn)
    params$tmax <- get_tmax(popn, params)

    message("Tmax = ", str_flatten_comma(params$tmax))
} else if (params$sim_new_data == "bici") {
    params$tmax <- c(t1 = 200, t2 = 200)
} else {
    # If Tinf, Tinc etc. are not in popn, create missing column(s) and set to NA
    missing_timings <- setdiff(params$timings, names(popn))
    if (length(missing_timings) > 0) {
        popn[, (missing_timings) := NA_real_]
        setcolorder(popn, params$timings, after = "donor")
    }
    params$tmax <- get_tmax(popn, params)
    
    # Correct for donors not properly infected
    popn <- fix_fb_data(popn, params)
    
    popn[Tdeath >= params$tmax[trial], Tdeath := NA]
}


### Plot the epidemic ----

plt <- plot_model(popn, params)



## Generate config files ----

{
    # Create missing directories
    params[str_ends(names(params), "_dir")] |>
        as.character() |>
        discard(dir.exists) |>
        walk(~ message("- mkdir ", .x)) |>
        walk(dir.create, recursive = TRUE)

    # Clean up old config files and generate fresh one
    cleanup_bici_files(params)
    bici_txt <- generate_bici_script(popn, params)
}


## Run BICI ----

{
    cmd <- with(params, str_glue(
        if (algorithm == "pas")
            "mpirun -n {nchains} --output :raw --oversubscribe " else "",
        "../BICI/bici-{platform} {data_dir}/{name}.bici {bici_cmd}",
        platform = Sys.info()[["sysname"]]
    ))
    message(str_glue("Running:\n$ {cmd}"))
    
    tic()
    out <- system(cmd)
    time_taken <- toc()
    
    if (out != 0) {
        stop("BICI failed to finish")
    }
}


## Retrieve results and save ----

{
    message("Retrieving results ...")
    
    name <- params$name
    dataset <- params$dataset
    output_dir <- params$output_dir
    results_dir <- params$results_dir
    bici_cmd <- params$bici_cmd
    
    if (bici_cmd == "inf") {
        # pe_name <- str_glue("{output_dir}/posterior.csv")
        # parameter_estimates <- fread(pe_name)
        parameter_estimates <- rebuild_bici_posteriors(dataset, name)
        
        ebvs_name     <- str_glue("{output_dir}/ebvs.csv")
        estimated_BVs <- if (file.exists(ebvs_name)) fread(ebvs_name)
        
        pa_name   <- str_glue("{output_dir}/pred_accs.csv")
        pred_accs <- if (file.exists(pa_name)) fread(pa_name)
        
        ranks <- if (!is.null(estimated_BVs)) get_ranks(popn, estimated_BVs, params)
        
        message("Parameter estimates:")
        msg_pars(parameter_estimates)
        
        results_pars <- c("params", "popn", "time_taken", "time_start", "time_end",
                          "parameter_estimates", "estimated_BVs","ranks", "pred_accs")
    } else {
        results_pars <- c("params", "popn", "time_taken", "time_start", "time_end")
    }
    
    time_end <- now()
    
    # Filter for results that we have and save
    results_pars |>
        keep(exists) |>
        mget() |>
        saveRDS(file = str_glue("{results_dir}/{name}.rds"))
    
    # Generate etc_inf.rds summary file
    flatten_bici_states(dataset, name, bici_cmd)
}

