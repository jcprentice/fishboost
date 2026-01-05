patch_in_traits <- function(popn, params) {
    message("Patching popn with traits ...")
    {
        setup           <- params$setup
        dataset         <- params$patch_dataset %||% ""
        name            <- params$patch_name
        model_traits    <- params$model_traits
        use_weight      <- params$use_weight
        patch_state     <- params$patch_state
        patch_type      <- params$patch_type
        traits_source   <- params$traits_source
        msgs            <- params$msgs
    }
    
    if (traits_source != "posterior" || dataset == "") {
        if (msgs) message(" - not requested")
        return(popn)
    }
    
    {
        base_dir    <- str_glue("datasets/{dataset}")
        data_dir    <- str_glue("{base_dir}/data")
        results_dir <- str_glue("{base_dir}/results")
        output_dir  <- str_glue("{data_dir}/{name}-out")
        states_dir  <- str_glue("{output_dir}/states")
        inf_out_dir <- str_glue("{output_dir}/output-inf")
    }
    
    # If we're not using mean, then need to sample from the files available, so
    # check if we even have any first
    if (patch_type %notin% c("mean", "median")) {
        f <- str_glue("{output_dir}/etc_inf.rds")

        # If we don't find anything, revert to patching with mean
        message(" - no ETC file, reverting to posterior mean")
        if (!file.exists(f)) patch_type <- "mean"
    }
    
    if (patch_type %in% c("mean", "median")) {
        f <- str_glue("{output_dir}/etc_inf.rds")
        if (file.exists(f)) {
            tmp <- readRDS(f)$popn # FIXME
            cols <- names(tmp) |> str_subset("^(sus|inf|lat|det|tol)")
            ebvs <- tmp[, map(.SD, mean), id, .SDcols = cols]
        } else {
            stop(str_glue(" - '{f}' not found!"))
        }
        
        if (msgs) message(str_glue(" - making popn by copying EBVs from '{dataset}/{name}' ..."))
        
    } else if (is.numeric(patch_state)) {
        f <- str_glue("{output_dir}/etc_inf.rds")
        if (file.exists(f)) {
            ebvs <- readRDS(f)$ies[state == patch_state, map(.SD, mean), id, .SDcols = -1]
        } else {
            stop(str_glue(" - '{f}' not found!"))
        }
        
        if (msgs) message(str_glue(" - making popn by copying EBVs from '{dataset}/{name}'\n",
                                   " - using state '{patch_state}' ..."))
    } else if (!is.null(params$trace_row)) {
        trace_row <- params$trace_row
        f <- str_glue("{output_dir}/extended_trace_combine.tsv")
        if (!file.exists(f)) {
            f <- str_remove(f, "extended_")
            if (!file.exists(f)) {
                stop("No trace file")
            }
        }
        if (msgs) message(str_glue(" - reading trace file '{f}'"))
        tc <- fread(f)[trace_row]
        tc <- tc[, .SD, .SDcols = str_subset(names(tc), "^\\d")] |>
            melt(measure.vars = measure(id, ie, ae, sep = "_"))
        ebvs <- tc[, .(id = as.integer(id),
                       ie = str_c(ie, "_", ae),
                       value)] |>
            dcast(id ~ ie)
        setnames(ebvs, c("id", str_c(names(ebvs)[2:7], " EBV")))
    } else {
        stop("No files to work with!")
    }
    
    if (nrow(popn) != nrow(ebvs)) {
        message(str_glue(" - inconsistent sizes for EBVs ({x}) and pedigree ({y}), generating new data",
                               x = nrow(popn), y = nrow(ebvs)))
        popn2 <- make_traits_from_pedigree(popn, params)
        return(popn2)
    }
    
    # Add missing traits (set to 0) and calculate phenotype
    walk(model_traits, \(x) {
        tg <- str_glue("{x}_g")
        te <- str_glue("{x}_e")
        if (tg %notin% names(ebvs)) {
            ebvs[, c(tg, te) := 0]
        }
        ebvs[, (x) := exp(get(tg) + get(te))]
    })
    
    # Reorder GT, EV, PT
    setcolorder(ebvs, c("id",
                        str_c(model_traits, "_g"),
                        str_c(model_traits, "_e"),
                        model_traits))
    
    # Now copy over to ebvs (I'm sure there must be an easier way to do this, but I
    # don't know what it is)
    popn2 <- copy(popn)
    walk(names(ebvs), ~ popn2[, (.x) := ebvs[, get(.x)]])
    
    popn2
}
