{
    source("libraries.R")
    source("source_files.R")
    source("get_Tinfs.R")
}

# Generate data sets ----
generate_km_data_sire <- function(dataset = "fb-final",
                             scen = 1,
                             opts = list(n_plots = 30,
                                         use_means = FALSE)) {
    
    # dataset <- "fb-final"; scen <- 1;
    # opts = list(n_plots = 3, = FALSE)
    # use_bici = TRUE
    n_plots      <- opts$n_plots
    use_means    <- opts$use_means
    
    message(str_glue("Generating KM data for '{dataset}/s{scen}-1'"))
    
    # Get parameters from target data set (and discard the bits we don't need)
    res_dir <- str_glue("datasets/{dataset}/results")
    files <- list.files(res_dir, pattern = str_glue("scen-{scen}-1")) |>
        str_sort(numeric = TRUE)
    
    if (is_empty(files)) stop(str_glue("No files for '{dataset}/s{scen}-1'"))
    
    fb_data <- map(files, \(f) {
        # f <- files[[1]]
        message(str_glue("Working with '{f}' ..."))
        res <- readRDS(str_glue("{res_dir}/{f}"))
        params <- res$params
        
        # We only need a subset of the fb_data
        popn <- res$popn
        popn[, `:=`(RP = Tdeath - Tsym, src = "fb")]
        popn[Tdeath == c(104, 160)[trial], Tdeath := NA]
        # popn <- fix_fb_data(popn)
        
        if ("parasites" %notin% names(popn)) {
            # parasites might be missing if we simulated
            popn[, parasites := !is.na(Tinf) | !is.na(Tsym) | !is.na(Tdeath)]
        }
        popn[donor == 1 & parasites == FALSE & is.na(Tdeath),
             `:=`(donor = 0, Tinf = NA)]
        popn[id == 1805, donor := 1]
        
        if ("GE" %notin% names(popn)) {
            GEs <- res$parameter_estimates[str_starts(parameter, "Group effect"), .(group = .I, GE = exp(mean))]
            popn <- merge(popn, GEs, by = "group", all = TRUE)
        }
        
        setcolorder(popn, "group", after = "trial")
        setcolorder(popn, "GE", after = "donor")
        
        list(popn = popn,
             params = params)
    })
    
    
    popns <- map(fb_data, \(fb) {
        # fb <- fb_data[[1]]
        params <- fb$params
        
        # Now override parameters
        params$sim_new_data    <- "bici"
        params$traits_source   <- "posterior"
        params$patch_dataset   <- dataset
        params$patch_name      <- str_glue("scen-{scen}-1")
        params$patch_type      <- if (use_means) "mean" else "sampled"
        params$patch_state     <- TRUE # str_starts(dataset, "fb")
        params$msgs <- FALSE
        
        # Check if weight1 and weight2 missing, if so add them
        if ("weight1" %notin% rownames(params$fe_vals)) {
            params$fe_vals <- with(params, rbind(fe_vals, fe_vals[rep("weight", 2),]))
            rownames(params$fe_vals)[5:6] <- str_c("weight", 1:2)
        }
        
        params_base <- copy(params)
        
        map(seq_len(n_plots), \(i) {
            # i <- 1L
            # This updates params with a different row from the trace file each time
            params <- params_base |>
                patch_params() |>
                # set_ge_opts() |>
                apply_links()
            
            
            ## Generate pedigree and popn ----
            # popn <- make_pedigree(params) |>
            #     set_groups(params) |>
            #     set_traits(params) |>
            #     apply_fixed_effects(params)
            popn <- fb$popn |>
                set_traits(params) |>
                apply_fixed_effects(params)
            
            # Simulate new data ----
            message(str_glue("{i} / {n_plots}"))
            
            popn <- simulate_epidemic(popn, params)
            popn[sdp == "progeny", .(sire, trial, donor, group,
                                     Tinf, Tsym, Tdeath,
                                     RP = Tdeath - Tsym, src = "sim")]
        })
    }) |>
        unlist(recursive = FALSE)
        
    message("Simulations complete, adding FB data")
    message()
    
    # Attach the actual data, labelling as "fb"
    fb <- map(fb_data, \(x) x$popn[sdp == "progeny", mget(names(popns[[1]]))])
    # NB: wrap fb_data in list()
    data <- rbindlist(c(popns, fb), fill = TRUE, idcol = "id")
    
    # Combine parts
    list(data = data,
         params = fb_data[[1]]$params,
         opts = c(dataset = dataset,
                  scen = scen,
                  opts))
}

