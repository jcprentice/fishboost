# if we want the same data as example.xml (handy for exploring the structure)

# sire_ex <- read_xml("example.xml")
# sire_ex_list <- as_list(sire_ex)

generate_sire22_xml <- function(data, params, GRM = NULL) {
    {
        name            <- params$name
        nsires          <- params$nsires
        nprogeny        <- params$nprogeny
        nparents        <- params$nparents
        ntotal          <- params$ntotal
        ngroups         <- params$ngroups
        model_type      <- params$model_type
        all_traitnames  <- params$all_traitnames
        traitnames      <- params$traitnames
        eta_type        <- params$eta_type
        rho_type        <- params$rho_type
        gamma_type      <- params$gamma_type
        group_effect    <- params$group_effect
        priors          <- params$priors
        time_step       <- params$time_step
        algorithm       <- params$algorithm
        anneal          <- params$anneal
        anneal_power    <- params$anneal_power
        nsample_per_gen <- params$nsample_per_gen
        nsample         <- params$nsample
        burnin          <- params$burnin
        thin            <- params$thin
        phi             <- params$phi
        data_dir        <- params$data_dir
        trial_fe        <- params$trial_fe
        donor_fe        <- params$donor_fe
        txd_fe          <- params$txd_fe
        link_traits     <- params$link_traits
        link_trial      <- params$link_trial
        link_donor      <- params$link_donor
        link_txd        <- params$link_txd
        link_shapes     <- params$link_shapes
        setup           <- params$setup
        use_fb_data     <- params$use_fb_data
        pass_events     <- params$pass_events
        tmax            <- params$tmax
        censor          <- params$censor
    }
    
    message(glue("Generating SIRE 2.2 XML file {data_dir}/{name}.xml ..."))
    
    # Begin XML file ----
    x <- xml_new_root("SIRE", version = "2.2")
    
    ## MCMC options ----
    xml_add_child(x, xml_comment("MCMC options"))
    
    if (algorithm == "pas") {
        xml_add_child(x, "pas", output_dir = glue("{data_dir}/{name}_out"),
                      nsample_per_gen = nsample_per_gen, phi_final = phi,
                      nsample = nsample, burnin = burnin, thin = thin)
    } else {
        xml_add_child(x, "mcmc", output_dir = glue("{data_dir}/{name}_out"),
                      nsample_per_gen = nsample_per_gen, phi_final = phi,
                      nsample = nsample, burnin = burnin, thin = thin,
                      anneal = anneal, anneal_power = anneal_power)
    }
    
    
    ## Model details ----
    
    # remove "_res" from model_type if it's there
    model_type <- sub("_res", "", model_type)
    
    ### Compartments ----
    
    # save pointer to compartment X tag as `cx` for adding fixed effects later
    
    xml_add_child(x, xml_comment("Model compartments"))
    if (grepl("S", model_type)) cs <- xml_add_child(x, "comp", name = "S")
    if (grepl("E", model_type)) ce <- xml_add_child(x, "comp", name = "E")
    if (grepl("I", model_type)) ci <- xml_add_child(x, "comp", name = "I", relative_infectivity = 1)
    if (grepl("D", model_type)) cd <- xml_add_child(x, "comp", name = "D", relative_infectivity = 1)
    if (grepl("R", model_type)) cr <- xml_add_child(x, "comp", name = "R")
    
    
    
    ### Transitions ----
    
    # link_shapes = "lll" means use "lat_shape" for det and rec
    shapes <- list(l = "eta_shape", d = "rho_shape", r = "gamma_shape")
    
    shape <- list(lat = shapes[[substr(link_shapes, 1, 1)]],
                  det = shapes[[substr(link_shapes, 2, 2)]],
                  rec = shapes[[substr(link_shapes, 3, 3)]])
    
    # save pointer to transition X tag as `tx` for adding fixed effects
    
    xml_add_child(x, xml_comment("Model transitions"))
    
    if (model_type == "SIR") {
        ts <- xml_add_child(x, "trans", from = "S", to = "I", type = "infection",beta = "beta", inf_model = "frequency dependent")
        tr <- xml_add_child(x, "trans", from = "I", to = "R", type = gamma_type, mean = "recovery_period", shape = shape$rec)
        if (time_step == 0) {
            xml_set_attr(ts, "data_column", "Tsym")
            xml_set_attr(tr, "data_column", "Trec")
        }
    } else if (model_type == "SEIR") {
        ts <- xml_add_child(x, "trans", from = "S", to = "E", type = "infection", beta = "beta", inf_model = "frequency dependent")
        tl <- xml_add_child(x, "trans", from = "E", to = "I", type = eta_type, mean = "latent_period", shape = shape$lat)
        tr <- xml_add_child(x, "trans", from = "I", to = "R", type = gamma_type, mean = "recovery_period", shape = shape$rec)
        if (time_step == 0) {
            if (pass_events == 3) xml_set_attr(ts, "data_column", "Tinf")
            xml_set_attr(tl, "data_column", "Tsym")
            xml_set_attr(tr, "data_column", "Trec")
        }
    } else if (model_type == "SIDR") {
        ts <- xml_add_child(x, "trans", from = "S", to = "I", type = "infection", beta = "beta", inf_model = "frequency dependent")
        td <- xml_add_child(x, "trans", from = "I", to = "D", type = rho_type, mean = "detection_period", shape = shape$det)
        tr <- xml_add_child(x, "trans", from = "D", to = "R", type = gamma_type, mean = "recovery_period", shape = shape$rec)
        if (time_step == 0) {
            if (pass_events >= 3) xml_set_attr(ts, "data_column", "Tinf")
            xml_set_attr(td, "data_column", "Tsym")
            xml_set_attr(tr, "data_column", "Trec")
        }
    } else if (model_type == "SEIDR") {
        ts <- xml_add_child(x, "trans", from = "S", to = "E", type = "infection", beta = "beta", inf_model = "frequency dependent")
        tl <- xml_add_child(x, "trans", from = "E", to = "I", type = eta_type, mean = "latent_period", shape = shape$lat)
        td <- xml_add_child(x, "trans", from = "I", to = "D", type = rho_type, mean = "detection_period", shape = shape$det)
        tr <- xml_add_child(x, "trans", from = "D", to = "R", type = gamma_type, mean = "recovery_period", shape = shape$rec)
        if (time_step == 0) {
            if (pass_events >= 4) xml_set_attr(ts, "data_column", "Tinf")
            if (pass_events >= 3) xml_set_attr(tl, "data_column", "Tinc")
            xml_set_attr(td, "data_column", "Tsym")
            xml_set_attr(tr, "data_column", "Trec")
        }
    }

    # Remove shape parameters when type is "exp"
    
    if (exists("tl") && eta_type == "exp") {
        xml_set_attr(tl, "shape", NULL)
        priors[parameter == "eta_shape", use := FALSE]
    }
    if (exists("td") && rho_type == "exp") {
        xml_set_attr(td, "shape", NULL)
        priors[parameter == "rho_shape", use := FALSE]
    }
    if (exists("tr") && gamma_type == "exp") {
        xml_set_attr(tr, "shape", NULL)
        priors[parameter == "gamma_shape", use := FALSE]
    }
    
    
    ### Add individual effects ----
    
    if ("susceptibility" %in% traitnames) {
        xml_set_attr(ts, "individual_effect", "s_a,s_e")
    }
    
    xml_add_child(x, xml_comment("Add FEs or IEs to infectivity here"))
    
    inf_node <- xml_add_child(x, "infectivity")
    if ("infectivity" %in% traitnames) {
        xml_set_attr(inf_node, "individual_effect", "i_a,i_e")
    }
    
    if ("recoverability" %in% traitnames) {
        xml_set_attr(tr, "individual_effect", "r_a,r_e")
        
        # link latency and detectability
        if ("latency" %in% traitnames && substr(link_traits, 2, 2) == "r") {
            xml_set_attr(tl, "individual_effect", "r_a,r_e")
        }
        
        if ("detectability" %in% traitnames && substr(link_traits, 4, 4) == "r") {
            xml_set_attr(td, "individual_effect", "r_a,r_e")
        }
    }
    
    
    
    
    ### Add fixed effects ----
    
    # Want to create vectors c(s="trial_s", ...) where the values may be linked
    all_traits1 <- substr(all_traitnames, 1, 1)
    
    trial_names <- paste0("trial_", all_traits1)
    zt <- match(strsplit(link_trial, "")[[1]], all_traits1)
    trial_names <- trial_names[zt]
    names(trial_names) <- all_traits1
    
    donor_names <- paste0("donor_", all_traits1)
    zd <- match(strsplit(link_donor, "")[[1]], all_traits1)
    donor_names <- donor_names[zd]
    names(donor_names) <- all_traits1
    
    txd_names <- paste0("txd_", all_traits1)
    zd <- match(strsplit(link_txd, "")[[1]], all_traits1)
    txd_names <- txd_names[zd]
    names(txd_names) <- all_traits1
    
    
    # Build a string for fixed_effects like "trial_i,donor_i,txd_i", and account
    # for linking of FEs
    fe_str <- function(x) {
        paste0(c(if (grepl(x, trial_fe)) trial_names[x] else NULL,
                 if (grepl(x, donor_fe)) donor_names[x] else NULL,
                 if (grepl(x, txd_fe))   txd_names[x]   else NULL),
               collapse = ",")
    }
    
    fex <- fe_str("s")
    if (fex != "") xml_set_attr(ts, "fixed_effect", fex)
    
    fex <- fe_str("l")
    if (fex != "") xml_set_attr(tl, "fixed_effect", fex)
    
    fex <- fe_str("i")
    if (fex != "") xml_set_attr(inf_node, "fixed_effect", fex)
    
    fex <- fe_str("d")
    if (fex != "") xml_set_attr(td, "fixed_effect", fex)
    
    fex <- fe_str("r")
    if (fex != "") xml_set_attr(tr, "fixed_effect", fex)
    
    
    ### Group effect ----
    if (group_effect < 0) {
        # useful if group_effect has been set outside of make_parameters()
        priors[parameter == "sigma", use := FALSE]
    } else {
        xml_add_child(x, xml_comment("Group effect"))
        xml_add_child(x, "group_effect", sigma = "sigma")
    }
    
    
    ## Covariance ----
    
    sir <- c("s", "i", "r")
    used_traits <- sir %in% substr(traitnames, 1, 1)
    names(used_traits) <- sir
    n_used <- sum(used_traits)
    use_cov <- any(used_traits)
    
    # Note that we have two child nodes both called "covariance", so need to
    # store the pointer, otherwise xml_add_child(x, "covariance") fails when we
    # try to add to the second one.
    
    if (use_cov) {
        # generate a string like "s_a,i_a" depending on which traits are used
        ie_a_str <- paste0(sir[used_traits], "_a", collapse = ",")
        ie_e_str <- gsub("a", "e", ie_a_str)
        
        cov_G_mat <- matrix(paste0("cov_G_", c("ss", "ii", "rr")[used_traits]), ncol = 1)
        cov_E_mat <- sub("G", "E", cov_G_mat)
        
        cor_G_mat <- diag(1, n_used, n_used)
        cors <- c(if (all(used_traits[c("s", "i")])) "r_G_si",
                  if (all(used_traits[c("s", "r")])) "r_G_sr",
                  if (all(used_traits[c("i", "r")])) "r_G_ir")
        cor_G_mat[upper.tri(cor_G_mat)] <- cors
        cor_G_mat[lower.tri(cor_G_mat)] <- t(cor_G_mat)[lower.tri(cor_G_mat)]
        cor_E_mat <- sub("G", "E", cor_G_mat)
        
        
        # Genetic Covariance
        xml_add_child(x, xml_comment("Genetic covariance between different individual effects"))
        xcg <- xml_add_child(x, "covariance",
                             individual_effect = ie_a_str,
                             relationship_matrix = "A")
        xml_add_child(xcg, "variance", table_to_tsv_string(cov_G_mat))
        xml_add_child(xcg, "correlation", table_to_tsv_string(cor_G_mat))
        
        # Environmental Covariance
        xml_add_child(x, xml_comment("Environmental covariance between different individual effects"))
        xce <- xml_add_child(x, "covariance",
                             individual_effect = ie_e_str,
                             relationship_matrix = "I")
        xml_add_child(xce, "variance", table_to_tsv_string(cov_E_mat))
        xml_add_child(xce, "correlation", table_to_tsv_string(cor_E_mat))
    }
    
    
    ## Observation and Inference periods ----
    xml_add_child(x, xml_comment("Inference and observation periods"))
    
    trials <- data[sdp == "progeny", unique(trial)]
    ntrials <- max(length(trials), 1L)
    for (i in seq_len(ntrials)) {
        groups <- data[trial == trials[i], mixedsort(unique(group))]
        trial_tmax <- if (use_fb_data | censor < 1) as.character(tmax[i]) else "infinity"
        
        for (tag in c("inference", "observation")) {
            xml_add_child(x, tag, group = paste(groups, collapse = ","), tmin = 0, tmax = trial_tmax)
        }
    }
    
    
    
    ## Priors ----
    xml_add_child(x, xml_comment("Model priors"))
    for (i in seq_len(nrow(priors))) {
        ppi <- as.list(priors[i, ])
        
        # skip over unused parameters
        if (ppi$use == FALSE) next
        
        xp <- xml_add_child(x, "prior")
        if (ppi$type == "Flat") {
            xml_set_attrs(xp, c(parameter = ppi$parameter, type = ppi$type, val1 = ppi$val1, val2 = ppi$val2))
        } else if (ppi$type == "Fixed") {
            xml_set_attrs(xp, c(parameter = ppi$parameter, type = ppi$type, val = ppi$val1))
        }
        # xml_add_child(x, "prior", parameter = ppi$parameter,
        #               type = "Flat", val1 = ppi$val1, val2 = ppi$val2)
    }
    
    
    ## Add the data table ----
    
    xml_add_child(x, xml_comment("Information about individuals"))
    
    dt <- glue("{data_dir}/{name}-data.tsv")
    write.table(data, file = dt, sep = "\t", quote = FALSE,
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
    if (!use_fb_data) {
        # for (x in c("s", "i", "r")) {
        #     if (x %in% substr(traitnames, 1, 1)) {
        #         xa = paste(x, "_a")
        #         xml_set_attr(xdt, xa, xa)
        #     }
        # }
        if (used_traits["s"]) xml_set_attr(xdt, "s_a", "s_a")
        if (used_traits["i"]) xml_set_attr(xdt, "i_a", "i_a")
        if (used_traits["r"]) xml_set_attr(xdt, "r_a", "r_a")
    }
    
    # add <datatable ... trial_r="9", donor_i="10", txd_i="11"> etc.
    
    
    for (str in unique(trial_names[strsplit(trial_fe, "")[[1]]])) {
        xml_set_attr(xdt, str, "trial_fe")
    }
    
    for (str in unique(donor_names[strsplit(donor_fe, "")[[1]]])) {
        xml_set_attr(xdt, str, "donor_fe")
    }
    
    for (str in unique(txd_names[strsplit(txd_fe, "")[[1]]])) {
        xml_set_attr(xdt, str, "txd_fe")
    }
    
    # for (i in strsplit(trial_fe, "")[[1]])
    # xml_set_attr(xdt, paste0("trial_", i), which(names(data) == "trial_fe"))
    # for (i in strsplit(donor_fe, "")[[1]])
    # xml_set_attr(xdt, paste0("donor_", i), which(names(data) == "donor_fe"))
    # for (i in strsplit(txd_fe, "")[[1]])
    # xml_set_attr(xdt, paste0("txd_", i), which(names(data) == "txd_fe"))
    
    
    ## Relationship matrix ----
    
    if (FALSE && !is.null(GRM)) {
        # Clip out the dams
        nondams <- c(seq.int(nsires), seq.int(nparents + 1L, ntotal))
        GRM2 <- GRM[nondams, nondams]
        
        # Convert the GRM into a sparse matrix, use `summary` to extract values,
        # then convert that into a data table
        
        GRM_dt <- as.data.table(summary(Matrix(GRM2, sparse = TRUE)))
        
        # Use the DT to make 0-based indices
        GRM_dt[, `:=`(i = i - 1L, j = j - 1L)]
        
        # Account for summary only giving UT matrix, and SIRE wanting access to
        # both the LT and UT parts
        GRM_dt <- rbind(GRM_dt, GRM_dt[i != j, .(i = j, j = i, x)])
        setkey(GRM_dt, j, i)
        
        
        # Write the DT as a TSV file
        A_file <- glue("{data_dir}/{name}-A.tsv")
        write.table(GRM_dt, file = A_file, sep = "\t", quote = FALSE,
                    row.names = FALSE, col.names = FALSE)

        xml_add_child(x, xml_comment("Relationship matrix (A or GRM)"))
        xml_add_child(x, "A_nonzero", file = A_file, name = "A")
    }
    
    ## Prediction Accuracies ----
    
    xml_add_child(x, xml_comment("Prediction Accuracies for sires"))
    xml_add_child(x, "prediction_accuracy", name = "sire",
                  ind = paste0(data[sdp == "sire", id], collapse = ","))
    
    xml_add_child(x, xml_comment("Prediction Accuracies for progeny"))
    xml_add_child(x, "prediction_accuracy", name = "progeny",
                  ind = paste0(data[sdp == "progeny", id], collapse = ","))
    
    
    # Write the XML file ----
    xml_file <- glue("{data_dir}/{name}.xml")
    write_xml(x, xml_file) #options = c("as_xml", "format") #?
    
    
    # Fix XML formatting ----
    
    # Seriously, R does a terrible job here, this makes the output much more readable.
    
    cmd <- paste0(
        # On Mac need to use GNU sed (gsed) instead of BSD sed
        if (Sys.info()["sysname"] == "Darwin") "g" else NULL,
        # sed: remove tabs at start; insert new lines before comments
        "sed -i 's/^\\t*//; s/^ *//; s/\\(<!--\\)/\\n\\1/' ",
        xml_file
    )
    system(cmd)
    
    return(x)
}
