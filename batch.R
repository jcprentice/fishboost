# Script for generating data and/or running inference in BICI on a server

## Load libraries and source files ----
suppressPackageStartupMessages({
    source("libraries.R")
    source("source_files.R")
})

time_start <- now()
message(str_glue("Started at {round(time_start)}"))

# set.seed(0)

## Set parameters ----

cmd_args <- commandArgs(trailingOnly = TRUE)
run_from_script <- length(cmd_args) > 0
if (run_from_script) {
    pname <- cmd_args[[1]]
    row_no <- as.integer(cmd_args[[2]])
} else {
    pname <- "fb-qtest"
    row_no <- 14L
}

{
    # Grab row_no from protocol file and common options (if they exist) and
    message(str_glue("Running: '{pname} / s{row_no}'"))
    protocol <- with(readRDS(str_glue("param_sets/{pname}.rds")),
                     safe_merge(protocol[row_no],
                                if (exists("common")) common))
    
    # protocol$patch_dataset <- "fb-test"
    
    # Create params with protocol
    params <- with(protocol, make_parameters(
        model_type = model_type,
        dataset = dataset,
        name = name,
        use_traits = use_traits,
        vars = vars,
        cors = cors,
        setup = setup,
        group_effect = group_effect,
        trial_fe = trial_fe,
        donor_fe = donor_fe,
        txd_fe = txd_fe,
        weight_fe = weight_fe
    ))
    
    message("Copying parameters ...")
    
    # This is necessary to ensure that `expand_priors` comes before any priors,
    # since the priors need to override it.
    protocol <- reorder_protocol(protocol)
    
    # Patch params from protocol file
    walk(names(protocol), \(param) {
        value <- protocol[[param]]
        
        # Skip NAs completely
        if (is.na(value)) return()
        
        if (is.character(value)) {
            message(str_glue("- {param} = '{value}'"))
        } else {
            message(str_glue("- {param} = {value}"))
        }
        
        # Skip any values already provided to make_pars()
        if (param %in% c("model_type", "dataset", "name", "setup",
                         "use_traits", "vars", "cors", "group_effect",
                         "trial_fe", "donor_fe", "txd_fe", "weight_fe")) {
            return()
        }
        
        # Handle special cases first
        if (str_detect(param, "periods")) {
            # set all of latent_period and latent_period_*,*
            sp <- str_split_i(param, "_", 1)
            params[[str_glue("{sp}_period")]] <<- value
            params$priors[str_detect(parameter, sp),
                          `:=`(true_val = value,
                               val2 = pmax(val2, value))]
            
        } else if (str_starts(param, "prior")) {
            # "prior__abc_xyz__val1" -> "abc_xyz", "val1"
            p1 <- str_split_1(param, "__")[-1]
            params$priors[parameter == p1[[1]], (p1[[2]]) := value]
            # params$priors[str_detect(parameter, p1[[1]]), (p1[[2]]) := value]
            
        } else if (param == "group_effect") {
            params[[param]] <<- value
            params$priors[parameter == "sigma", `:=`(
                val2 = max(2 * value, 1e-3),
                true_val = max(value, 0))]
            
        } else if (str_starts(param, "fe_")) {
            # fe_trial_inf -> fe_vals["trial", "inf"]
            p1 <- str_split_1(param, "_")[-1]
            td <- p1[[1]]; tr <- p1[[2]]; tr1 <- str_sub(tr, 1, 1)
            params$fe_vals[[td, tr]] <<- value
            params$priors[parameter == str_c(td, "_", tr1),
                          `:=`(val1 = min(value - 2, val1),
                               val2 = max(value + 2, val2),
                               true_val = value)]
            
        } else if (str_ends(param, "_fe")) {
            # These should already have been passed to make_params(), and don't
            # want to overwrite "" with "none" or NA
            return()
            
        } else if (param == "weight_is_nested") {
            params$weight_is_nested <<- if (params$ntrials == 1) {
                FALSE
            } else {
                value
            }
            
        } else if (param == "expand_priors") {
            params$priors[str_detect(parameter, "trial|donor|txd|weight"),
                          `:=`(val1 = -value, val2 = value)]
            
        } else if (param == "pass_events") {
            params$pass_events <<- if (value == "all") {
                params$timings
            } else {
                str_split_1(value, ",")
            }
            
        } else if (param == "ge_opts") {
            params$ge_opts <<- str_split_1(value, ",")
            
        } else if (param == "fix_donors") {
            params$fix_donors <<- str_split_1(value, ",")
            
        } else if (param == "t_demote") {
            vals <- as.numeric(str_split_1(as.character(value), ","))
            params$t_demote <<- if (length(vals) == 1) rep(vals, 2) else vals
            
        } else if (param == "skip_patches") {
            # By default patch_params will try to patch everything. This makes
            # it skip certain sections such as covariances.
            params$skip_patches <<- str_split_1(value, ",")
            
        } else if (param == "h2") {
            params$cov_G <<- params$cov_G * value
            params$cov_E <<- params$cov_E * (1 - value)
            # FIXME: do something about Sigma_X
            
        } else {
            # Any other parameter that doesn't require special treatment
            params[param] <<- value
        }
    })
    
    # Tidy up numbers (only in balanced / non-Fishboost cases)
    if (params$sim_new_data != "no" && !str_starts(params$setup, "fb")) {
        params$ndams      <- with(params, nsires * dpsire)
        params$nparents   <- with(params, nsires + ndams)
        params$nprogeny   <- with(params, ndams * ppdam)
        params$ntotal     <- with(params, nparents + nprogeny)
        params$group_size <- with(params, nprogeny / ngroups)
    }
    
    # Overrule whatever is in the protocol and just take it from qsub
    num_slots <- as.integer(Sys.getenv("NUM_SLOTS"))
    if (!is.na(num_slots) && params$nchains != num_slots) {
        message("- Overriding nchains with $NUM_SLOTS = ", num_slots)
        params$nchains <- num_slots
    }
    
    if (!run_from_script) {
        # For testing purposes
        params$nchains <- 2
        params$nsample <- 1e3
    }
    
    if (params$bici_cmd %in% c("sim", "post-sim")) {
        params$nchains <- with(params, min(nchains, nreps))
    }
    
    # This can go awry
    if (is.na(params$group_effect)) {
        params$group_effect <- -1
    }
    
    # Set seed before any randomness is used
    tmp_seed <- protocol$seed %||% 0
    params$seed <- if (tmp_seed >= 0) tmp_seed else protocol$replicate
    set.seed(params$seed)
    
    params2 <- params |>
        patch_params() |>  # Patch params with posteriors from dataset / scenario
        # set_ge_opts() |>   # Genetic covariance options
        set_use_flags() |> # Ensure priors are correctly enabled
        apply_links()      # Fix any traits that need to be linked
    params <- params2
    
    # Testing, remove when checked
    if (params$single_prior == "uniform") {
        params$priors[type == "inverse", type := "uniform"]
    }
    
    # Tidy up LP, DP, RP, including rate, shape, and priors
    params <- tidy_up_periods(params, protocol)
    
    params$R0 <- with(params, r_beta * (detection_period + removal_period))
    
    # params$DEBUG <- TRUE
    
    # Quick check of how things are looking
    summarise_params(params)
}

## Generate pedigree and popn ----

if (params$sim_new_data %in% c("r", "bici")) {
    popn <- make_pedigree(params) |>
        set_groups(params) |>
        set_traits(params) |>
        set_weights(params) |>
        apply_fixed_effects(params)
} else if (params$sim_new_data == "no") {
    # Load popn, pedigree, and GRM
    popn <- readRDS(str_glue("fb_data/{params$setup}.rds")) |>
        fix_fb_data(params)
    
    params$tmax <- c(t1 = 104, t2 = 160)
    
    # FB only has the last 2 events
    params$pass_events <- c("Tsym", "Tdeath")
} else {
    rf <- with(params, str_glue("datasets/{patch_dataset}/data/",
                                "{patch_name}-out/{sim_new_data}.rds"))
    popn <- readRDS(rf)$popn[state == max(params$patch_state, 1L)]
    popn[, state := NULL]
}


## Run epidemic ----

if (params$sim_new_data == "r") {
    popn <- simulate_epidemic(popn, params)
    popn[sdp == "progeny", parasites := !is.na(Tinf)]
    params$estimated_R0 <- get_R0(popn)
    params$tmax <- get_tmax(popn, params)
} else if (params$sim_new_data == "bici") {
    params$tmax <- c(t1 = 200, t2 = 200)
} else {
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

message(str_glue("censor = {x}\nTmax = {y}",
                 x = params$censor,
                 y = str_flatten_comma(params$tmax)))


### Plot the epidemic ----

# plt <- plot_model(popn, params)


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
        "../BICI/bici-{platform} {config}.bici {bici_cmd}",
        platform = Sys.info()[["sysname"]]
    ))
    message(str_glue("Running:\n$ {cmd}"))
    
    nattempts <- if (run_from_script) 4 else 1
    for (attempt in seq_len(nattempts)) {
        tic()
        out <- system(cmd)
        time_taken <- toc()
        
        if (out == 0) {
            message("BICI ran successfully")
            break
        } else if (attempt < nattempts) {
            message(str_glue("- attempt {attempt}/{nattempts} failed, trying again"))
            inc_seed(params)
        } else {
            stop("BICI failed to finish")
        }
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
        # parameter_estimates <- if (file.exists(pe_name)) fread(pe_name)
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
    
    # Generate etc_inf.rds or etc_sim.rds summary file
    flatten_bici_states(dataset, name, bici_cmd)
}

message("Finished!")
