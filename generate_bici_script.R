# if we want the same data as example.xml (handy for exploring the structure)

# Takes a list and converts it into a BICI string, and wraps strings in "", e.g.
# comp name="S" x=1 y=2 z="fnoof"

generate_bici_script <- function(popn, params) {
    with(params, {
        # attach(params)
        message(str_glue("Generating BICI script file '{data_dir}/{name}.bici' ..."))
        
        if (DEBUG) {
            # This turns off overwriting files while debugging
            message("DEBUG is on! No files will be written.")
            fwrite <- function(...) {
                message(str_glue("fwrite to '{..2}'"))
                print(..1)
            }
            writeLines <- function(...) {
                message(str_glue("writeLines to '{..2}'"))
                print(..1)
            }
        }
        
        out_dir <- str_glue("{data_dir}/{name}-out")
        if (!dir.exists(out_dir)) {
            message("- mkdir ", out_dir)
            dir.create(out_dir)
        }
        
        # Begin BICI script ----
        x <- list()
        
        # Data directory ----

        x$section_description <- "Data directory"

        x$data_dir <- list(node = "data-dir", folder = str_glue("{name}-out"))
        
        x$description <- list(node = "description", text = "description.txt")
        
        writeLines(c(str_glue("Dataset: {dataset} / {name}"),
                     str_glue("File written: {round_date(now())}"),
                     description),
                   str_glue("{out_dir}/description.txt"))
        
        
        # Main setup ----
        
        x$section_setup <- "Setup"
        
        ## Inference ----
        
        x$comment_inf = "Inference options"
        
        seed <- as.integer(seed) %% 1e4L
        timestep <- max(ceiling(max(tmax) / 10) / 1e3, time_step_bici)
            
        x$inference <- list(node = "inference",
                            diagnostics = "on",
                            start = 0,
                            end = max(tmax),
                            timestep = timestep,
                            seed = seed)
        
        # algorithm specific inference opts
        alg <- if (algorithm == "pas") {
            list(algorithm = "PAS-MCMC",
                 update = as.integer(nsample),
                 npart = as.integer(nchains),
                 # "npart-per-core" = 1L,
                 "gen-percent" = 2L,
                 "param-output" = as.integer(min(nsample, thinto)),
                 "state-output" = max(as.integer(min(nsample, sample_states)), 10L),
                 "burnin-percent" = 100 * burnprop)
            
        } else if (algorithm == "da") {
            list(algorithm = "DA-MCMC",
                 update = as.integer(nsample),
                 nchain = as.integer(nchains),
                 "chain-per-core" = as.integer(nchains),
                 "param-output" = as.integer(min(nsample, thinto)),
                 "state-output" = max(sample_states, 10L),
                 "burnin-percent" = 100 * burnprop,
                 anneal = "none")
            
        } else if (algorithm == "abc") {
            list(algorithm = "ABC",
                 sample = 1000L,
                 "acc-frac" = 0.1)
            
        } else if (algorithm == "abc-smc") {
            list(algorithm = "ABC-SMC",
                 sample = 1000L,
                 "acc-frac" = 0.5,
                 gen = 5,
                 "kernel-size" = 0.5)
            
        } else if (algorithm == "pmcmc") {
            list(algorithm = "PMCMC")
            
        } else if (algorithm == "mfa") {
            list(algorithm = "MFA")
        }
        
        x$inference <- c(x$inference, alg)
        
        
        ## Simulation ----
        
        x$comment_sim <- "Simulation options"
        
        tmax10 <- ceiling(max(tmax) / 10) * 10
        
        x$simulation <- list(node = "simulation",
                             start = 0,
                             end = tmax10,
                             number = nreps,
                             timestep = timestep,
                             seed = seed)
       
        ## Post-Sim ----
        
        x$comment_postsim <- "Posterior simulation options"
        
        x$postsim <- list(node = "post-sim",
                          start = 0,
                          end = tmax10,
                          number = nreps,
                          seed = seed)
        
        
        # Model details ----
        
        x$section_model_species <- "Model species"
        
        x$species <- list(node = "species", name = "Fish", type = "individual")
        
        
        ## Compartments ----
        
        x$comment_model_compartments <- "Model compartments"
        
        ### Disease status ----
        
        x$class_ds <- list(node = "class", name = "DS", index = "a")
        
        x$comp_ds_S <- list(node = "comp", name = "S")
        x$comp_ds_E <- list(node = "comp", name = "E")
        x$comp_ds_I <- list(node = "comp", name = "I")
        x$comp_ds_D <- list(node = "comp", name = "D")
        x$comp_ds_R <- list(node = "comp", name = "R")
        
        # Remove nodes for models we don't need
        if (model_type %in% c("SIR", "SIDR")) x$comp_ds_E <- NULL
        if (model_type %in% c("SIR", "SEIR")) x$comp_ds_D <- NULL
        
        # The flashy way:
        # walk(model_type |> uniq_chars(),
        #      ~ {x[[str_glue("comp_ds_{.x}")]] <<- list(node = "comp", name = .x)})
        
        ### Transitions ----
        
        x$comment_transitions <- "Model transitions"
        
        # Function to build transition period and infection time distributions.
        # This includes place holders for individual, fixed, and group effects,
        # which can be string replaced later.
        
        # Transition period
        get_TP <- function(p = "l", dist_type = "exp") {
            m <- c(l = "LP", d = "DP", t = "RP")[[p]]
            # link_shapes = "lll" means use "eta^shape" for det and tol
            cv1 <- c(l = "LP^shape", d = "DP^shape", t = "RP^shape")
            ldt <- names(cv1)
            cv <- setNames(cv1[str_chars(link_shapes)], ldt)[[p]]
            
            switch(dist_type, 
                   "gamma" = str_glue("gamma(mean:FE_IE_{m}_b,c, cv:{cv}_b,c)"),
                   "exp"   = str_glue("exp(mean:FE_IE_{m}_b,c)"))
        }
        
        # Infection time
        if (!exists("inf_model")) inf_model <- 1
        priors[parameter == "infrat", use := inf_model == 4]
        
        get_IT <- function(id = "I|D") {
            # this is rate, because it's a poisson process, not a waiting time
            inf_eqn <- if (id == "I|D" || inf_model != 2) {
                switch(inf_model,
                       str_glue("{{{id},g; FEi_IEi_}}"),
                       str_glue("(0.1*{{I,g; FEi_IEi_}}+{{D,g; FEi_IEi_}})"),
                       str_glue("(0.1*{{{id},g,Don; FEi_IEi_}+{{{id},g,Rec; FEi_IEi_}})"),
                       str_glue("(infrat*{{{id},g,Don; FEi_IEi_}+{{{id},g,Rec; FEi_IEi_}})"))
            } else {
                if (inf_model == 2) {
                    message("Warning, incompatible `inf_model`, defaulting to `inf_model=1`")
                }
                str_glue("{{{id},g; FEi_IEi_}}/{group_size}")
            }
            str_glue("exp(rate:GE_FEs_IEs_beta_b{inf_eqn}/{group_size})")
        }
        
        if (model_type == "SIR") {
            x$trans_inf <- list(node = "trans", name = "S->I", value = get_IT("I"))
            x$trans_tol <- list(node = "trans", name = "I->R", value = get_TP("l", RP_dist))
        } else if (model_type == "SEIR") {
            x$trans_inf <- list(node = "trans", name = "S->E", value = get_IT("I"))
            x$trans_lat <- list(node = "trans", name = "E->I", value = get_TP("l", LP_dist))
            x$trans_tol <- list(node = "trans", name = "I->R", value = get_TP("t", RP_dist))
        } else if (model_type == "SIDR") {
            x$trans_inf <- list(node = "trans", name = "S->I", value = get_IT("I|D"))
            x$trans_det <- list(node = "trans", name = "I->D", value = get_TP("d", DP_dist))
            x$trans_tol <- list(node = "trans", name = "D->R", value = get_TP("t", RP_dist))
        } else if (model_type == "SEIDR") {
            x$trans_inf <- list(node = "trans", name = "S->E", value = get_IT("I|D"))
            x$trans_lat <- list(node = "trans", name = "E->I", value = get_TP("l", LP_dist))
            x$trans_det <- list(node = "trans", name = "I->D", value = get_TP("d", DP_dist))
            x$trans_tol <- list(node = "trans", name = "D->R", value = get_TP("t", RP_dist))
        }
        
        
        # Trial, Donor, Group ----
        
        ## Trial ----
        x$comment_trial_compartments <- "Trial"
        
        x$class_trial <- list(node = "class", name = "Trial", index = "b")
        
        trials <- popn[sdp == "progeny", unique(trial)]
        walk(trials, ~ {
            x[[str_glue("comp_trial_{.x}")]] <<-
                list(node = "comp", name = str_glue("Tr{.x}"))
        })
        
        
        ## Donor ----
        x$comment_donor_compartments <- "Donor"
        
        x$class_donor <- list(node = "class", name = "Donor", index = "c")
        x$comp_group_don <- list(node = "comp", name = "Don")
        x$comp_group_rec <- list(node = "comp", name = "Rec")
        
        
        ## Group ----
        x$comment_group_compartments <- "Group"
        
        x$class_group <- list(node = "class", name = "Group", index = "g")
        x$comp_group <- list(node = "comp-all", fix = "true", file = "comp-group.tsv")
        
        groups <- popn[!is.na(group), sort(unique(group))]
        N <- ceiling(sqrt(length(groups)))
        data.table(name = groups,
                   expand.grid(x = 10 * seq(0, N - 1),
                               y = 10 * seq(0, N - 1))[seq_along(groups), ]) |>
            fwrite(file = str_glue("{out_dir}/comp-group.tsv"),
                   sep = "\t", quote = TRUE)
        
        
        # Data ----
        
        x$comment_data <- "Data for inference and simulation"
        
        x$add_ind_inf <- list(node = "add-ind-inf", file = "add-ind.tsv")
        x$rem_ind_inf <- list(node = "remove-ind-inf", file = "remove-ind.tsv")
        
        popn[sdp == "progeny", .(ID = id,
                                 t = 0,
                                 Group = group,
                                 Trial = str_c("Tr", trial),
                                 Donor = fifelse(donor == 0, "Rec", "Don"),
                                 DS = compartments[1 + donor])] |>
            fwrite(file = str_glue("{out_dir}/add-ind.tsv"),
                   sep = "\t", quote = TRUE)
        
        popn[sdp == "progeny", .(ID = id,
                                 t = tmax[str_c("t", trial)])] |>
            fwrite(file = str_glue("{out_dir}/remove-ind.tsv"),
                   sep = "\t", quote = TRUE)
        
        if (sim_new_data != "bici") {
            if (popn_format == "intervals") {
                message("- Writing intervals")
                x$comp_data <- list(node = "comp-data",
                                    name = "DS observations",
                                    class = "DS",
                                    file = "comp-data-DS.tsv")
                
                discretise_time_bici(popn, params) |>
                    fwrite(file = str_glue("{out_dir}/comp-data-DS.tsv"),
                           sep = "\t", quote = TRUE)
                
            } else if (popn_format == "times") {
                message("- Writing times")
                
                # Censor the data
                data <- popn[, .SD, .SDcols = c("id", "trial", "donor", timings)]
                data[, (names(.SD)) := map(.SD, \(x) fifelse(x > tmax[trial], NA, x, NA)),
                     trial, .SDcols = timings]
                
                walk(seq_along(timings), \(i) {
                    ti <- timings[[i]]
                    if (ti %notin% pass_events) return()
                    
                    e1 <- compartments[[i]]
                    e2 <- compartments[[i + 1]]
                    file_str <- str_glue("trans-data-{e1}{e2}.tsv")
                    xtd <- list(node = "trans-data",
                                name = "Transition data",
                                trans = str_glue("{e1}->{e2}"),
                                obsrange = "all",
                                file = file_str)
                    
                    x[[str_glue("trans_data_{e1}{e2}")]] <<- xtd
                    
                    out <- if (ti == "Tinf" && ti %in% pass_events) {
                        data[donor == 0 & !is.na(Tinf), .(ID = id, t = Tinf)]
                    } else {
                        data[!is.na(get(ti)), .(ID = id, t = get(ti))]
                    }
                    fwrite(out, file = str_glue("{out_dir}/{file_str}"), sep = "\t", quote = TRUE)
                })
                
                rm(data)
            }
        }
        
        
        x$comment_simulation <- "Simulation data"
        
        x$add_ind_sim <- list(node = "add-ind-sim", file = "add-ind.tsv")
        x$rem_ind_sim <- list(node = "remove-ind-sim", file = "remove-ind.tsv")
        
        # x$comment_post_simulation <- "Post-simulation data"
        # 
        # x$add_ind_post_sim <- list(node = "add-ind-post-sim", file = "add-ind.tsv")
        # x$rem_ind_post_sim <- list(node = "remove-ind-post-sim", file = "remove-ind.tsv")
        
        
        x$section_effects <- "Individual / fixed / group effects"
        
        ## Individual Effects ----
        
        sildt <- c("s", "i", "l", "d", "t")
        
        
        x$comment_ies <- "Individual effects"
        
        # For visualising groups of individuals in BICI
        
        x$sire <- list(node = "ind-group-data", name = "sire", file = "sire.tsv")
        popn[sdp == "sire", .(ID = id)] |>
            fwrite(file = str_glue("{out_dir}/sire.tsv"), sep = "\t", quote = TRUE)
        
        x$dam <- list(node = "ind-group-data", name = "dam", file = "dam.tsv")
        popn[sdp == "dam", .(ID = id)] |>
            fwrite(file = str_glue("{out_dir}/dam.tsv"), sep = "\t", quote = TRUE)
        
        x$progeny <- list(node = "ind-group-data", name = "progeny", file = "progeny.tsv")
        popn[sdp == "progeny", .(ID = id)] |>
            fwrite(file = str_glue("{out_dir}/progeny.tsv"), sep = "\t", quote = TRUE)
        
        
        ## Additive genetic effects ----
        
        # IEs as "sg,ig,tg". We use this several times, so keep it.
        ies_str <- str_c(intersect(str_chars(use_traits),
                                   str_chars(link_traits)),
                         "g", collapse = ",")
    
        if (use_traits %notin% c("", "none")) {
            
            # BICI can only handle one of inverse or sparse, prefer inverse
            if (str_detect(use_grm, "H.*inv.*nz")) {
                use_grm <- str_remove(use_grm, "_nz")
            }
                
            x$comment_grm <- if (str_starts(use_grm, "A")) {
                "Using sparse A matrix"
            } else if (str_starts(use_grm, "H")) {
                str_glue("Using{s}{i} {H} matrix",
                         s = if (str_detect(use_grm, "nz")) " sparse" else "",
                         i = if (str_detect(use_grm, "inv")) " inverse" else "",
                         H = str_split_i(use_grm, "_", 1))
            } else {
                "Using pedigree"
            }
            
            x$ie_g <- list(node = "ind-effect",
                           name = "gen",
                           ie = ies_str)
            
            message(str_glue("- using GRM = '{use_grm}'"))
            
            if (use_grm %in% c("", "none", "pedigree")) {
                # Write the pedigree file
                xgrm <- list(pedigree = "pedigree.tsv")
                
                popn[, .(ID = as.character(id),
                         sire = fifelse(is.na(sire), ".", as.character(sire), "."),
                         dam  = fifelse(is.na(dam),  ".", as.character(dam),  "."))] |>
                    fwrite(file = str_glue("{out_dir}/pedigree.tsv"),
                           sep = "\t", quote = TRUE)
                
            } else if (str_starts(use_grm, "A")) {
                # Write the A_nz matrix regardless of what use_grm is
                GRM <- make_grm(popn, "A_nz")
                xgrm <- list("A-sparse" = "A-matrix.tsv",
                             "ind-list" = "ind-list.tsv")
                popn[, .(Individual = id)] |>
                    fwrite(str_glue("{out_dir}/ind-list.tsv"),
                           sep = "\t", quote = TRUE)
                
                if (min(GRM$i) == 1) GRM[, `:=`(i = i - 1, j = j - 1)]
                
                fwrite(GRM, file = str_glue("{out_dir}/A-matrix.tsv"), sep = "\t")
                
            } else {
                # Write the pre-computed H matrix
                
                H_mat <- str_glue("{use_grm}_{str}.tsv",
                                  str = str_remove(setup, "fb_"))
                
                xgrm <- if (str_detect(use_grm, "nz")) {
                    popn[, .(Individual = id)] |>
                        fwrite(str_glue("{out_dir}/ind-list.tsv"),
                               sep = "\t", quote = TRUE)
                    list("A-sparse" = H_mat,
                         "ind-list" = "ind-list.tsv")
                } else if (str_detect(use_grm, "inv")) {
                    list(Ainv = H_mat)
                } else {
                    list(A = H_mat)
                }
                
                # Force relative soft link
                system(str_glue("{g}ln -sfr fb_data/{H_mat} {out_dir}/{H_mat}",
                                # Use GNU ln not BSD ln
                                g = if (Sys.info()[["sysname"]] == "Darwin")
                                    "g" else ""))
            }
            
            x$ie_g <- c(x$ie_g, xgrm)
            
            # Environmental effects
            x$ie_e <- list(node = "ind-effect",
                           name = "env",
                           ie = str_replace_all(ies_str, "a|g", "e"))
        }
        
        ## True BVs ----
        if (bici_cmd == "inf" && str_detect(dataset, "sim")) {
            traits_to_fit <- model_traits[intersect(str_chars(use_traits),
                                                    str_chars(link_traits))]
            ie_names <- expand.grid(traits_to_fit, c("g", "e")) |>
                apply(1, str_flatten, "_")
            ie_vars <- expand.grid(str_1st(traits_to_fit),
                                   c("g", "e")) |> apply(1, str_flatten)
            
            if (FALSE) {
                # All in one file `ie-data.tsv`
                walk(seq_along(ie_names), \(i) {
                    ie_name <- ie_names[[i]]
                    ie_var <- ie_vars[[i]]
                    
                    # ind_effect-data name="sus_g" ie="sg" file="ie-data.tsv" cols="ID,sg"
                    x[[str_glue("ie_{ie_var}")]] <<- list(
                        node = "ind-effect-data",
                        name = ie_name,
                        ie = ie_var,
                        file = "ie-data.tsv",
                        cols = str_glue("ID,{ie_var}"))
                })
                
                popn[, .SD, .SDcols = c("id", ie_names)] |>
                    setnames(c("id", ie_names), c("ID", ie_vars)) |>
                    fwrite(file = str_glue("{out_dir}/ie-data.tsv"),
                           sep = "\t")
            } else {
                # Split between files `ie-sg-data.tsv`
                walk(seq_along(ie_names), \(i) {
                    ie_name <- ie_names[[i]]
                    ie_var <- ie_vars[[i]]
                    fstr <- str_glue("ie-{ie_var}-data.tsv")
                    
                    # ind-effect-data name="sus_g" ie="sg" file="ie-sg-data.tsv"
                    x[[str_glue("ie_{ie_var}")]] <<- list(
                        node = "ind-effect-data",
                        name = ie_name,
                        ie = ie_var,
                        file = fstr)
                    
                    popn[, .(ID = id,
                             Value = if (ie_name %in% names(popn)) get(ie_name) else 0)] |>
                        fwrite(file = str_glue("{out_dir}/{fstr}"),
                               sep = "\t")
                })
            }
        }
        
        
        ## Fixed effects ----
        
        x$comment_fe <- "Fixed effects"
        
        # Get the appropriate linked effect for a given trait
        get_LE <- function(links, t1) setNames(str_chars(links), sildt)[[t1]]
        
        get_FEs <- function(trait) {
            t1 <- str_1st(trait)
            fe <- if (str_detect(weight_fe, t1)) {
                fex <- get_LE(link_weight, t1)
                if (weight_is_nested && ntrials > 1) {
                    str_glue("<weight1{fex}><weight2{fex}>")
                } else {
                    str_glue("<weight{fex}>")
                }
            } else ""
            str_flatten(fe)
        }
        
        get_IEs <- function(trait) {
            t1 <- str_1st(trait)
            if (str_detect(use_traits, t1)) {
                iex <- get_LE(link_traits, t1)
                str_glue("[{iex}g][{iex}e]")
            } else {
                ""
            }
        }
        
        if ("trans_inf" %in% names(x)) {
            x$trans_inf$value <- x$trans_inf$value |>
                str_replace_all(c("GE_" = if (group_effect >= 0) "G_g*" else "",
                                  "FEs_" = get_FEs("sus"), "IEs_" = get_IEs("sus"),
                                  "FEi_" = get_FEs("inf"),
                                  "IEi_" = get_IEs("inf")))
        }
        
        walk(c("lat", "det", "tol"), \(trait) {
            tx <- str_glue("trans_{trait}")
            if (tx %in% names(x)) {
                x[[tx]]$value <<- x[[tx]]$value |>
                    str_replace_all(c("FE_" = get_FEs(trait),
                                      "IE_" = get_IEs(trait)))
            }
        })
        
        # If no FEs or Ies on infection, then remove the "; "
        x$trans_inf$value <- str_replace_all(x$trans_inf$value, "; \\}", "}")
        
        ### Weight ----
        
        wt_fes <- priors[str_starts(parameter, "weight") & use, parameter]
        
        
        walk(wt_fes, ~ {
            wt <- str_split_i(.x, "_", 1)
            x[[str_glue("fe_{.x}")]] <<- list(
                node = "fixed-effect",
                name = str_remove(.x, "_"),
                X = str_c(str_glue("fe-{wt}.tsv")))
        })
        
        if ("weight" %in% names(popn)) {
            rec_fn <- if (use_weight == "log") log_recentre else recentre
            if (weight_is_nested) {
                popn[, `:=`(weight1_fe = 0, weight2_fe = 0)]
                popn[trial == 1, weight1_fe := rec_fn(weight)]
                popn[trial == 2, weight2_fe := rec_fn(weight)]
            } else {
                popn[, weight_fe := 0]
                popn[sdp == "progeny", weight_fe := rec_fn(weight)]
            }
        }
        
        
        if (weight_is_nested) {
            pop2 <- popn[sdp == "progeny", .(ID = id, value = log_recentre(weight)), trial]
            walk(1:2, \(i) {
                pop2[, .(ID, value = fifelse(trial == i, value, 0))] |>
                    fwrite(str_glue("{out_dir}/fe-weight{i}.tsv"),
                           sep = "\t", quote = TRUE)
            })
            rm(pop2)
        } else {
            popn[sdp == "progeny", .(ID = id, value = log_recentre(weight))] |>
                fwrite(str_glue("{out_dir}/fe-weight.tsv"),
                       sep = "\t", quote = TRUE)
        }
        
        
        
        ## Group effect ----
        if (group_effect < 0) {
            # useful if group_effect has been set outside of make_parameters()
            priors[parameter == "sigma", use := FALSE]
        } else {
            # x$section_group_effect <- "Group effect"
            # x$group_effect <- list(sigma = "sigma")
        }
        
        
        # Priors ----
        
        x$section_priors <- "Model priors"
        
        period_to_pp <- function(period = "latent") {
            period |> str_1st() |> str_to_upper() |> str_c("P")
        }
        
        # Check for nested weights
        if (weight_is_nested && ntrials > 1) {
            priors[str_starts(parameter, "weight_"), use := FALSE]
        } else {
            priors[str_starts(parameter, "weight[12]_"), use := FALSE]
        }
        
        # Check use for gamma distributed LPs
        priors[parameter == "LP_shape", use := LP_dist == "gamma"]
        priors[parameter == "DP_shape", use := DP_dist == "gamma"]
        priors[parameter == "RP_shape", use := RP_dist == "gamma"]
        
        # BICI treats FEs as classes, remove them as FEs, and add appropriate
        # files like 'value-beta.tsv'. We don't have to use them if they aren't
        # needed by the FE.
        
        {
            priors[str_starts(parameter, "trial|donor|txd"), use := FALSE]
            
            trials <- popn[sdp == "progeny", str_c("Tr", unique(trial))]
            
            # beta
            bp <- priors[parameter %in% str_c("beta_", trials),
                         .(Value = true_val,
                           Prior = str_glue("{single_prior}(0.01,{x})", x = val2))]
            
            data.table(b = trials, Value = bp$Value) |>
                fwrite(str_glue("{out_dir}/value-beta.tsv"),
                       sep = "\t", quote = TRUE)
            
            data.table(b = trials, Prior = bp$Prior) |>
                fwrite(str_glue("{out_dir}/prior-beta.tsv"),
                       sep = "\t", quote = TRUE)
            
            
            tab <- expand.grid(c = c("Don", "Rec"), b = trials) |> rev()
            t_bc <- with(tab, str_c("_", b, ",", c))
            
            # Shape parameters
            walk(c("LP", "DP", "RP"), \(sp) {
                if (get(str_glue("{sp}_dist")) == "exp") return()
                
                Value <- priors[parameter == str_glue("{sp}_shape"), true_val]
                data.table(tab, Value) |>
                    fwrite(str_glue("{out_dir}/value-{sp}^shape.tsv"),
                           sep = "\t", quote = TRUE)
            })
            
            
            # Event periods
            walk(c("latent", "detection", "removal"), \(period) {
                pp <- period_to_pp(period)
                
                prior <- priors[parameter %in% str_c(period, "_period", t_bc),
                                .(Value = true_val,
                                  Prior = fcase(
                                      type == "constant",
                                      str_glue("fix({x})", x = true_val),
                                      type == "uniform",
                                      str_glue("uniform({v1},{v2})",v1 = val1, v2 = val2),
                                      type == "inverse",
                                      str_glue("inverse({v1},{v2})", v1 = max(1, val1), v2 = val2)
                                  )), .I]
                
                
                data.table(tab, Value = prior$Value) |>
                    fwrite(str_glue("{out_dir}/value-{pp}.tsv"),
                           sep = "\t", quote = TRUE)
                
                data.table(tab, Prior = prior$Prior) |>
                    fwrite(str_glue("{out_dir}/prior-{pp}.tsv"),
                           sep = "\t", quote = TRUE)
            })
        }
        
        # Check for fixed periods
        lp_types <- with(priors[str_ends(parameter, "period"), .(parameter, type)],
                         setNames(type, parameter)) |> as.list()
        
        if (use_traits %notin% c("", "none")) {
            val1 <- params$priors[str_detect(parameter, "cov_"), min(val1)] |>
                max(if (cov_prior == "uniform") 0.01 else 0.5)
            
            val2 <- params$priors[str_detect(parameter, "cov_"), max(val2)]
            
            if (!exists("cov_prior")) cov_prior <- "uniform"
            
            x$prior_cov_G <- list(
                node = "param",
                name = "\\Omega^gen_z,z'",
                value = "value-cov-gen.tsv",
                prior = str_glue("mvn-{cov_prior}({val1},{val2})")
            )
            
            # Possibly already defined, maybe not
            traits_to_fit <- model_traits[intersect(str_chars(use_traits),
                                                    str_chars(link_traits))]
            
            XG <- Sigma_G[traits_to_fit, traits_to_fit] |> round(5)
            dimnames(XG) <- list(NULL, str_c(str_1st(traits_to_fit), "g"))
            XG[lower.tri(XG)] <- "."
            
            as.data.table(XG) |>
                fwrite(str_glue("{out_dir}/value-cov-gen.tsv"),
                       sep = "\t", quote = TRUE)
            
            x$prior_cov_E <- list(
                node = "param",
                name = "\\Omega^env_z,z'",
                value = "value-cov-env.tsv",
                prior = str_glue("mvn-{cov_prior}({val1},{val2})")
            )
                
            XE <- Sigma_E[traits_to_fit, traits_to_fit] |> round(5)
            dimnames(XE) <- list(NULL, str_c(str_1st(traits_to_fit), "e"))
            XE[lower.tri(XE)] <- "."
            
            as.data.table(XE) |>
                fwrite(str_glue("{out_dir}/value-cov-env.tsv"),
                       sep = "\t", quote = TRUE)
        }
        
        walk(seq_len(nrow(priors)), \(i) {
            ppi <- as.list(priors[i])
            
            # skip over unused parameters
            if (ppi$use == FALSE) return()
            
            # Ensure numbers like 1e-3 appear as 0.001
            ppi <- map(ppi, ~ format(.x, scientific = FALSE))
            
            pname <- ppi$parameter
            node_str <- str_glue("prior_{pname}")
            
            if (pname == "beta") {
                pname <- "beta_b"
                ppi$true_val <- "value-beta.tsv"
            } else if (str_ends(pname, "period")) {
                # "latent_period" -> "LP_b,c"
                pname <- period_to_pp(pname)
                ppi$true_val <- str_glue("value-{pname}.tsv")
                pname <- pname |> str_c("_b,c")
            } else if (str_ends(pname, "shape")) {
                # "LP_shape" -> "LP^shape_b,c"
                pname <- pname |> str_replace("_shape", "^shape")
                ppi$val1 <- 0.45
                ppi$val2 <- 1.4
                ppi$true_val <- str_glue("value-{pname}.tsv")
                pname <- pname |> str_c("_b,c")
            } else if (str_detect(pname, "sigma")) {
                # Note: group effect gets 2 entries
                ppi$val1 <- 1e-2
                x$prior_ge <<- list(node = "param", name = "G_g",
                                    dist = "log-normal(1,sigma)")
                ppi$type <- "inverse"
            } else if (str_starts(pname, "weight")) {
                # "weight1-s" -> "\mu^weight1s"
                pname <- pname |> str_remove("_") |> str_c("\\mu^", x = _)
            } else if (str_detect(pname, "_[GE]_")) {
                return()
            }
            # } else if (str_starts(pname, "cov_[GE]_")) {
            #     # "cov_G_ss" -> "Ω^sa,sa"
            #     tmp <- str_split_1(pname, "_")
            #     ge <- str_to_lower(tmp[[2]])
            #     tr <- tmp[[3]] |> str_1st()
            #     pname <- str_glue("\\Omega^{tr}{ge},{tr}{ge}")
            # } else if (str_starts(pname, "r_[GE]_")) {
            #     # "r_G_si" -> "ω^sg,ig"
            #     tmp <- str_split_1(pname, "_")
            #     ge <- str_to_lower(tmp[[2]])
            #     tr <- str_chars(tmp[[3]])
            #     pname <- str_glue("\\omega^{tr[[1]]}{ge},{tr[[2]]}{ge}")
            # }
            
            x_ns <- list(node = "param", name = pname, value = ppi$true_val)
            
            if (ppi$parameter == "beta") {
                x_ns[["prior-split"]] <- "prior-beta.tsv"
            } else if (str_ends(ppi$parameter, "period")) {
                x_ns[["prior-split"]] <- str_glue("prior-{ppar}.tsv",
                                                  ppar = period_to_pp(ppi$parameter))
            } else {
                x_ns$prior <- switch(
                    ppi$type,
                    "constant" = str_glue("fix({ppi$true_val})"),
                    "uniform"  = str_glue("uniform({ppi$val1},{ppi$val2})"),
                    "inverse"  = str_glue("inverse({v1},{ppi$val2})",
                                          v1 = max(ppi$val1, 0.01))
                )
            }
            
            x[[node_str]] <<- x_ns
        })
        
        if (FALSE) {
            x[str_detect(names(x), "prior")] |> str()
        }
        
        
        # Write the BICI files ----
        
        # Pad section headers with empty lines on either side for readability,
        # and fix comments
        sections <- str_subset(names(x), "^section")
        x[sections] <- str_c("\n\n# -------- ", str_to_upper(x[sections]), " --------\n")
        # Remove empty lines from start
        x[[1]] <- str_remove(x[[1]], "\\n+")
        
        comments <- str_subset(names(x), "^comment")
        x[comments] <- str_c("\n# ", x[comments])
        
        bici_str <- function(node) {
            # If the node is just a character, not a list, leave as is
            if (is.character(node))
                return(as.character(node))
            
            # Strings have to be wrapped in ""
            wrap_str <- function(x) if (is.character(x)) str_c('"', x, '"') else as.character(x)
            # wrap_str <- function(x) str_c('"', as.character(x), '"')
            
            map2_chr(names(node), node, \(x, y) {
                if (x == "node") y else str_c(x, "=", wrap_str(y))
                # if (is.null(x) || x == "") y else str_c(x, " = ", wrap_str(y))
            }) |>
                str_c(collapse = " ")
        }
        
        x1 <- map_chr(x, bici_str) |> unname()
        
        writeLines(x1, str_glue("{data_dir}/{name}.bici"))
        
        # detach(params)
        invisible(x)
    })
}
