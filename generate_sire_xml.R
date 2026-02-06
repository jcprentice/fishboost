library(xml2)

# if we want the same data as example.xml (handy for exploring the structure)

# sire_ex <- read_xml("example.xml")
# sire_ex_list <- as_list(sire_ex)

generate_sire_xml <- function(popn, params) {
    with(params, {
        # attach(params)
        message(str_glue("Generating SIRE 2.2 XML file {data_dir}/{name}.xml ..."))
        
        # Begin XML file ----
        x <- xml_new_root("SIRE", version = "2.2")
        
        ## MCMC options ----
        xml_add_child(x, xml_comment("MCMC options"))
        
        xalg <- xml_add_child(x, algorithm,
                              output_dir = str_glue("{data_dir}/{name}-out"),
                              phi_final = phi,
                              nsample = as.integer(nsample),
                              # burnin = as.integer(burnin),
                              burnin = as.integer(ceiling(nsample * burnprop)),
                              # thin = as.integer(thin),
                              thin = as.integer(ceiling(nsample / thinto)),
                              # nsample_per_gen = as.integer(nsample_per_gen),
                              nsample_per_gen = as.integer(max(2.5e-3 * nsample, 1)),
                              sample_states = sample_states, ie_output = ie_output)
        if (algorithm == "mcmc") {
            xml_set_attrs(xalg, c(anneal = anneal, anneal_power = anneal_power))
        }
        
        
        ## Model details ----
        
        # remove "_res" from model_type if it's there
        model_type <- str_remove(model_type, "_res")
        
        ### Compartments ----
        
        # save pointer to compartment X tag as `cx` for adding fixed effects later
        
        xml_add_child(x, xml_comment("Model compartments"))
        
        cs <- xml_add_child(x, "comp", name = "S")
        ce <- xml_add_child(x, "comp", name = "E")
        ci <- xml_add_child(x, "comp", name = "I", relative_infectivity = 1)
        cd <- xml_add_child(x, "comp", name = "D", relative_infectivity = 1)
        cr <- xml_add_child(x, "comp", name = "R")
        
        if (model_type %in% c("SIR", "SIDR")) xml_remove(ce)
        if (model_type %in% c("SIR", "SEIR")) xml_remove(cd)
        
        
        
        ### Transitions ----
        
        # link_shapes = "lll" means use "lat_shape" for det and tol
        shapes <- list(l = "LP_shape", d = "DP_shape", t = "RP_shape")
        
        shape <- list(lat = shapes[[str_sub(link_shapes, 1, 1)]],
                      det = shapes[[str_sub(link_shapes, 2, 2)]],
                      tol = shapes[[str_sub(link_shapes, 3, 3)]])
        
        # save pointer to transition X tag as `tx` for adding fixed effects
        
        xml_add_child(x, xml_comment("Model transitions"))
        
        if (model_type == "SIR") {
            ts <- xml_add_child(x, "trans", from = "S", to = "I", type = "infection", beta = "beta", inf_model = "frequency dependent")
            tr <- xml_add_child(x, "trans", from = "I", to = "R", type = RP_dist, mean = "removal_period", shape = shape$tol)
            if (time_step == 0) {
                if ("Tsym" %in% pass_events) xml_set_attr(ts, "data_column", "Tsym")
                if ("Tdeath" %in% pass_events) xml_set_attr(tr, "data_column", "Tdeath")
            }
        } else if (model_type == "SEIR") {
            ts <- xml_add_child(x, "trans", from = "S", to = "E", type = "infection", beta = "beta", inf_model = "frequency dependent")
            tl <- xml_add_child(x, "trans", from = "E", to = "I", type = LP_dist, mean = "latent_period", shape = shape$lat)
            tr <- xml_add_child(x, "trans", from = "I", to = "R", type = RP_dist, mean = "removal_period", shape = shape$tol)
            if (time_step == 0) {
                if ("Tinf" %in% pass_events) xml_set_attr(ts, "data_column", "Tinf")
                if ("Tsym" %in% pass_events) xml_set_attr(tl, "data_column", "Tsym")
                if ("Tdeath" %in% pass_events) xml_set_attr(tr, "data_column", "Tdeath")
            }
        } else if (model_type == "SIDR") {
            ts <- xml_add_child(x, "trans", from = "S", to = "I", type = "infection", beta = "beta", inf_model = "frequency dependent")
            td <- xml_add_child(x, "trans", from = "I", to = "D", type = DP_dist, mean = "detection_period", shape = shape$det)
            tr <- xml_add_child(x, "trans", from = "D", to = "R", type = RP_dist, mean = "removal_period", shape = shape$tol)
            if (time_step == 0) {
                if ("Tinf" %in% pass_events) xml_set_attr(ts, "data_column", "Tinf")
                if ("Tsym" %in% pass_events) xml_set_attr(td, "data_column", "Tsym")
                if ("Tdeath" %in% pass_events) xml_set_attr(tr, "data_column", "Tdeath")
            }
        } else if (model_type == "SEIDR") {
            ts <- xml_add_child(x, "trans", from = "S", to = "E", type = "infection", beta = "beta", inf_model = "frequency dependent")
            tl <- xml_add_child(x, "trans", from = "E", to = "I", type = LP_dist, mean = "latent_period", shape = shape$lat)
            td <- xml_add_child(x, "trans", from = "I", to = "D", type = DP_dist, mean = "detection_period", shape = shape$det)
            tr <- xml_add_child(x, "trans", from = "D", to = "R", type = RP_dist, mean = "removal_period", shape = shape$tol)
            if (time_step == 0) {
                if ("Tinf" %in% pass_events) xml_set_attr(ts, "data_column", "Tinf")
                if ("Tinc" %in% pass_events) xml_set_attr(tl, "data_column", "Tinc")
                if ("Tsym" %in% pass_events) xml_set_attr(td, "data_column", "Tsym")
                if ("Tdeath" %in% pass_events) xml_set_attr(tr, "data_column", "Tdeath")
            }
        }
        
        # Remove shape parameters when distribution is "exp"
        
        if (exists("tl") && LP_dist == "exp") {
            xml_set_attr(tl, "shape", NULL)
            priors[parameter == "LP_shape", use := FALSE]
        }
        if (exists("td") && DP_dist == "exp") {
            xml_set_attr(td, "shape", NULL)
            priors[parameter == "DP_shape", use := FALSE]
        }
        if (exists("tr") && RP_dist == "exp") {
            xml_set_attr(tr, "shape", NULL)
            priors[parameter == "RP_shape", use := FALSE]
        }
        
        
        ### Add individual effects ----
        sildt <- c("s", "i", "l", "d", "t")
        traits1 <- str_sub(traitnames, 1, 1)
        
        # evs to remove
        no_evs <- if (any(str_starts(ge_opts, "no_ev"))) {
            ge_opts |> str_subset("no_ev") |> str_remove("no_ev_") |> str_chars()
        } else if ("gt_only" %in% ge_opts) {
            sildt
        }
        
        link_traits1 <- str_chars(link_traits)
        xaxe <- if ("pt_only" %in% ge_opts) {
            str_c(link_traits1, "_e") |> setNames(sildt)
        } else if ("gt_only" %in% ge_opts) {
            str_c(link_traits1, "_g") |> setNames(sildt)
        } else if (any(str_starts(ge_opts, "no_ev"))) {
            str_c(link_traits1, "_g,", link_traits1, "_e") |>
                str_remove(str_c(",", no_evs, "_e", collapse = "|")) |>
                setNames(sildt)
        } else {
            str_c(link_traits1, "_g,", link_traits1, "_e") |> setNames(sildt)
        }
        
        xml_add_child(x, xml_comment("Add FEs or IEs to infectivity here"))
        ti <- xml_add_child(x, "infectivity")
        
        if ("s" %in% traits1) xml_set_attr(ts, "individual_effect", xaxe[["s"]])
        if ("i" %in% traits1) xml_set_attr(ti, "individual_effect", xaxe[["i"]])
        if ("l" %in% traits1) xml_set_attr(tl, "individual_effect", xaxe[["l"]])
        if ("d" %in% traits1) xml_set_attr(td, "individual_effect", xaxe[["d"]])
        if ("t" %in% traits1) xml_set_attr(tr, "individual_effect", xaxe[["t"]])
        
        
        ### Add fixed effects ----
        
        # Want to create vectors c(s="trial_s", ...) where the values may be linked
        trial_names   <- setNames(str_c("trial_",   str_chars(link_trial)),  sildt)
        donor_names   <- setNames(str_c("donor_",   str_chars(link_donor)),  sildt)
        txd_names     <- setNames(str_c("txd_",     str_chars(link_txd)),    sildt)
        weight_names  <- setNames(str_c("weight_",  str_chars(link_weight)), sildt)
        weight1_names <- setNames(str_c("weight1_", str_chars(link_weight)), sildt)
        weight2_names <- setNames(str_c("weight2_", str_chars(link_weight)), sildt)
        
        
        # Build a string for fixed_effects like "trial_i,donor_i,txd_i", and account
        # for linking of FEs
        fe_str <- function(x) {
            str_flatten(c(if (str_detect(trial_fe,  x)) trial_names[[x]],
                          if (str_detect(donor_fe,  x)) donor_names[[x]],
                          if (str_detect(txd_fe,    x)) txd_names[[x]],
                          if (str_detect(weight_fe, x)) {
                              if (weight_is_nested) {
                                  c(weight1_names[[x]], weight2_names[[x]])
                              } else {
                                  weight_names[[x]]
                              }
                          }),
                        ",")
        }
        
        fex <- fe_str("s")
        if (exists("ts") && fex != "") xml_set_attr(ts, "fixed_effect", fex)
        
        fex <- fe_str("i")
        if (exists("ti") && fex != "") xml_set_attr(ti, "fixed_effect", fex)
        
        fex <- fe_str("l")
        if (exists("tl") && fex != "") xml_set_attr(tl, "fixed_effect", fex)
        
        fex <- fe_str("d")
        if (exists("td") && fex != "") xml_set_attr(td, "fixed_effect", fex)
        
        fex <- fe_str("t")
        if (exists("tr") && fex != "") xml_set_attr(tr, "fixed_effect", fex)
        
        
        ### Group effect ----
        if (group_effect < 0) {
            # useful if group_effect has been set outside of make_parameters()
            priors[parameter == "sigma", use := FALSE]
        } else {
            xml_add_child(x, xml_comment("Group effect"))
            xml_add_child(x, "group_effect", sigma = "sigma")
        }
        
        
        ## Covariance ----
        
        sit <- c("s", "i", "t")
        used_traits <- setNames(sit %in% traits1, sit)
        used_traits_e <- setNames(sit %in% traits1 & sit %notin% no_evs, sit)
        n_used <- sum(used_traits)
        n_used_e <- sum(used_traits_e)
        
        
        # Note that we have two child nodes both called "covariance", so need to
        # store the pointer, otherwise xml_add_child(x, "covariance") fails when we
        # try to add to the second one.
        
        if (n_used > 0) {
            # generate a string like "s_g,i_g" depending on which traits are used
            ie_g_str <- str_c(sit[used_traits], "_g", collapse = ",")
            ie_e_str <- str_c(sit[used_traits_e], "_e", collapse = ",")
            
            cov_G_mat <- matrix(str_c("cov_G_", c("ss", "ii", "tt")[used_traits]), ncol = 1)
            cov_E_mat <- matrix(str_c("cov_E_", c("ss", "ii", "tt")[used_traits_e]), ncol = 1)
            
            cor_G_mat <- diag(1, n_used, n_used)
            cors <- c(if (all(used_traits[c("s", "i")])) "r_G_si",
                      if (all(used_traits[c("s", "t")])) "r_G_st",
                      if (all(used_traits[c("i", "t")])) "r_G_it")
            cor_G_mat[upper.tri(cor_G_mat)] <- cors
            cor_G_mat[lower.tri(cor_G_mat)] <- t(cor_G_mat)[lower.tri(cor_G_mat)]
            
            cor_E_mat <- diag(1, n_used_e, n_used_e)
            cors <- c(if (all(used_traits_e[c("s", "i")])) "r_E_si",
                      if (all(used_traits_e[c("s", "t")])) "r_E_st",
                      if (all(used_traits_e[c("i", "t")])) "r_E_it")
            cor_E_mat[upper.tri(cor_E_mat)] <- cors
            cor_E_mat[lower.tri(cor_E_mat)] <- t(cor_E_mat)[lower.tri(cor_E_mat)]
            
            ### Genetic Covariance ----
            
            # Convert matrices to TSV strings
            cov_G_tsv <- table_to_tsv_string(cov_G_mat)
            cov_E_tsv <- table_to_tsv_string(cov_E_mat)
            cor_G_tsv <- table_to_tsv_string(cor_G_mat)
            cor_E_tsv <- table_to_tsv_string(cor_E_mat)
            
            if ("pt_only" %notin% ge_opts) {
                xml_add_child(x, xml_comment("Genetic covariance between different individual effects"))
                xcg <- xml_add_child(x, "covariance",
                                     individual_effect = ie_g_str,
                                     relationship_matrix = if (str_starts(use_grm, "H")) "H" else "A")
                xml_add_child(xcg, "variance", cov_G_tsv)
                xml_add_child(xcg, "correlation", cor_G_tsv)
            }
            
            ### Environmental Covariance ----
            if ("gt_only" %notin% ge_opts || n_used_e > 0) {
                xml_add_child(x, xml_comment("Environmental covariance between different individual effects"))
                xce <- xml_add_child(x, "covariance",
                                     individual_effect = ie_e_str,
                                     relationship_matrix = "I")
                xml_add_child(xce, "variance", cov_E_tsv)
                xml_add_child(xce, "correlation", cor_E_tsv)
            }
        }
        
        
        ## Observation and Inference periods ----
        xml_add_child(x, xml_comment("Inference and observation periods"))
        
        walk(popn[sdp == "progeny", unique(trial)], \(i) {
            groups <- popn[trial == i, group] |>
                unique() |>
                str_flatten(",")
            
            walk(c("inference", "observation"), \(tag) {
                xml_add_child(x, tag, group = groups,
                              tmin = 0, tmax = tmax[[str_c("t", i)]])
            })
        })
        
        
        ## Priors ----
        
        # Check for nested weights
        if (weight_is_nested) {
            priors[str_detect(parameter, "weight_"), use := FALSE]
        } else {
            priors[str_detect(parameter, "weight[12]_"), use := FALSE]
        }
        
        xml_add_child(x, xml_comment("Model priors"))
        
        walk(seq_len(nrow(priors)), \(i) {
            ppi <- as.list(priors[i, ])
            
            # skip over unused parameters
            if (ppi$use == FALSE)
                return()
            
            xp <- xml_add_child(x, "prior")
            if (ppi$type == "constant") {
                xml_set_attrs(xp, c(parameter = ppi$parameter, type = "Fixed", val = ppi$true_val))
            } else {
                xml_set_attrs(xp, c(parameter = ppi$parameter, type = ppi$type, val1 = ppi$val1, val2 = ppi$val2))
            }
        })
        
        
        ## Add the data table ----
        xml_add_child(x, xml_comment("Information about individuals"))
        
        dt <- str_glue("{data_dir}/{name}-data.tsv")
        if (msgs) message("- Writing data table to ", dt)
        
        prepare_data(popn, params) |>
            write.table(file = dt, sep = "\t", quote = FALSE,
                        row.names = FALSE, col.names = TRUE)
        
        xdt <- xml_add_child(x, "datatable", file = dt,
                             id = "id", sire = "sire", dam = "dam", group = "group",
                             initial_comp = "initial_comp", prediction_accuracy = "sdp")
        if (time_step > 0) {
            xml_set_attr(xdt, "comp_status", "comp_status")
        }
        
        # Write the data to a TSV file, then use sed to replace the string
        # "REPLACE_WITH_DATA". This is trickier, but *much* faster than trying
        # to generate a huge string in the XML stream.
        
        
        
        
        ## DT col names ----
        if (sim_new_data != "no") {
            if (used_traits["s"]) xml_set_attr(xdt, "s_g", "s_g")
            if (used_traits["i"]) xml_set_attr(xdt, "i_g", "i_g")
            if (used_traits["t"]) xml_set_attr(xdt, "t_g", "t_g")
        }
        
        walk(unique(trial_names[str_chars(trial_fe)]), \(str) {
            xml_set_attr(xdt, str, "trial_fe")
        })
        
        walk(unique(donor_names[str_chars(donor_fe)]), \(str) {
            xml_set_attr(xdt, str, "donor_fe")
        })
        
        walk(unique(txd_names[str_chars(txd_fe)]), \(str) {
            xml_set_attr(xdt, str, "txd_fe")
        })
        
        if (weight_is_nested) {
            walk(unique(weight1_names[str_chars(weight_fe)]), \(str) {
                xml_set_attr(xdt, str, "weight1_fe")
            })
            
            walk(unique(weight2_names[str_chars(weight_fe)]), \(str) {
                xml_set_attr(xdt, str, "weight2_fe")
            })
        } else {
            walk(unique(weight_names[str_chars(weight_fe)]), \(str) {
                xml_set_attr(xdt, str, "weight_fe")
            })
        }
        
        
        ## Genomic Relationship Matrix ----
        
        if (use_grm %in% c("", "none", "pedigree")) {
            # No GRM
        } else if (str_starts(use_grm, "A")) {
            GRM <- make_grm(popn, use_grm)
            
            A_file <- str_glue("{data_dir}/{name}-A.tsv")
            write.table(GRM, file = A_file, sep = "\t", quote = FALSE,
                        row.names = FALSE, col.names = FALSE)
            
            if (str_detect(use_grm, "inv")) {
                xml_add_child(x, xml_comment("Inverse pedigree matrix (A_inv)"))
                xml_add_child(x, "A_inv_nz", file = A_file, name = "A")
            } else {
                xml_add_child(x, xml_comment("Pedigree matrix (A)"))
                xml_add_child(x, "A_nz", file = A_file, name = "A")
            }
            # xml_add_child(x, "A_nonzero", file = A_file, name = "A")
            if (msgs) message("- Writing GRM to ", A_file)
            
        } else if (str_starts(use_grm, "H")) {
            
            if (str_starts(setup, "fb")) {
                # Just point to the pre-existing H matrix for FB
                str <- str_remove(setup, "fb_")
                sire_H_file <- str_glue("fb_data/{use_grm}_{str}.tsv")
                if (msgs) message(str_glue("- pointing {use_grm} to {sire_H_file}"))
            } else {
                # Write the passed GRM
                write.table(GRM, file = sire_H_file, sep = "\t", quote = FALSE,
                            row.names = FALSE, col.names = FALSE)
                sire_H_file <- str_glue("{data_dir}/{name}-H.tsv")
                if (msgs) message(str_glue("- Writing GRM to {sire_H_file}"))
            }
            
            if (str_detect(use_grm, "inv")) {
                xml_add_child(x, xml_comment("Inverse H matrix"))
                xml_add_child(x, "Ainv", file = sire_H_file, name = "H")
            } else {
                xml_add_child(x, xml_comment("H matrix"))
                xml_add_child(x, "A", file = sire_H_file, name = "H")
            }
        }
        
        
        ## Prediction Accuracies ----
        
        xml_add_child(x, xml_comment("Prediction Accuracies for sires"))
        xml_add_child(x, "prediction_accuracy", name = "sire",
                      ind = str_flatten(popn[sdp == "sire", id], ","))
        
        xml_add_child(x, xml_comment("Prediction Accuracies for dams"))
        xml_add_child(x, "prediction_accuracy", name = "dam",
                      ind = str_flatten(popn[sdp == "dam", id], ","))
        
        xml_add_child(x, xml_comment("Prediction Accuracies for progeny"))
        xml_add_child(x, "prediction_accuracy", name = "progeny",
                      ind = str_flatten(popn[sdp == "progeny", id], ","))
        
        
        # Write the XML file ----
        xml_file <- str_glue("{data_dir}/{name}.xml")
        write_xml(x, xml_file) #options = c("as_xml", "format") #?
        
        
        # Fix XML formatting ----
        
        # R does a terrible job here, this makes the output much more readable.
        # On Mac need to use GNU sed (gsed) instead of BSD sed.
        
        cmd <- str_c(
            if (Sys.info()[["sysname"]] == "Darwin") "g", # force gnu sed
            "sed -i - e'",          # sed, in place, run script
            "s/^\\s*//; ",          # remove white space from start of line
            "s/\\(<!--\\)/\\n\\1/", # add newline before comment
            "' ",                   # end script
            xml_file                # /path/to/file.xml
        )
        system(cmd)
        
        invisible(x)
    })
}
