# if we want the same data as example.xml (handy for exploring the structure)

# sire_ex <- read_xml("example.xml")
# sire_ex_list <- as_list(sire_ex)

generate_sire21_xml <- function(data, params, GRM) {
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
        group_effect    <- params$group_effect
        priors          <- params$priors
        time_step       <- params$time_step
        nsample         <- params$nsample
        burnin          <- params$burnin
        thin            <- params$thin
        anneal          <- params$anneal
        anneal_power    <- params$anneal_power
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
    }
    
    message(glue("Generating SIRE 2.1 XML file {data_dir}/{name}.xml ..."))
    
    # Begin XML file ----
    x <- xml_new_root("SIRE", version = "2.1")
    
    ## MCMC options ----
    xml_add_child(x, xml_comment("MCMC options"))
    xml_add_child(x, "mcmc", output_dir = glue("{data_dir}/{name}_out"),
                  nsample = nsample, burnin = burnin, thin = thin,
                  anneal = anneal, anneal_power = anneal_power)
    
    
    ## Model details ----
    
    # remove "_res" from model_type if it's there
    model_type <- sub("_res", "", model_type)
    
    ### Compartments ----
    
    # save pointer to compartment X tag as `cx` for adding fixed effects
    
    xml_add_child(x, xml_comment("Model compartments"))
    if (grepl("S", model_type)) cs <- xml_add_child(x, "comp", name = "S")
    if (grepl("E", model_type)) ce <- xml_add_child(x, "comp", name = "E")
    if (grepl("I", model_type)) ci <- xml_add_child(x, "comp", name = "I", infectivity = 1)
    if (grepl("D", model_type)) cd <- xml_add_child(x, "comp", name = "D", infectivity = 1)
    if (grepl("R", model_type)) cr <- xml_add_child(x, "comp", name = "R")
    
    
    ### Transitions ----
    
    # link_shapes = "lll" means use "lat_shape" for det and rec
    shapes <- list(l = "eta_shape", d = "rho_shape", r = "gamma_shape")
    shape <- list(lat = shapes[[substr(link_shapes, 1, 1)]],
                  det = shapes[[substr(link_shapes, 2, 2)]],
                  rec = shapes[[substr(link_shapes, 3, 3)]])
    
    # save pointer to transition X tag as `tx` for adding fixed effects
    
    xml_add_child(x, xml_comment("Model transitions"))

	message("time_step = ", time_step)
    
    switch(
        model_type,
        "SIR" = {
            ts <- xml_add_child(x, "trans", from = "S", to = "I", type = "infection", beta = "beta", inf_model = "frequency dependent")
            tr <- xml_add_child(x, "trans", from = "I", to = "R", type = "gamma", mean = "recovery_period", shape = shape$rec)
            if (time_step == 0) {
                xml_set_attr(ts, "data_column", 4)
                xml_set_attr(tr, "data_column", 5)
            }
        },
        "SEIR" = {
            ts <- xml_add_child(x, "trans", from = "S", to = "E", type = "infection", beta = "beta", inf_model = "frequency dependent")
            tl <- xml_add_child(x, "trans", from = "E", to = "I", type = "gamma", mean = "latent_period", shape = shape$lat)
            tr <- xml_add_child(x, "trans", from = "I", to = "R", type = "gamma", mean = "recovery_period", shape = shape$rec)
            if (time_step == 0) {
                pass_events <- clamp(pass_events, 2, 3)
                if (pass_events == 3) xml_set_attr(ts, "data_column", pass_events + 1)
                xml_set_attr(tl, "data_column", pass_events + 2)
                xml_set_attr(tr, "data_column", pass_events + 3)
            }
        },
        "SIDR" = {
            ts <- xml_add_child(x, "trans", from = "S", to = "I", type = "infection", beta = "beta", inf_model = "frequency dependent")
            td <- xml_add_child(x, "trans", from = "I", to = "D", type = "gamma", mean = "detection_period", shape = shape$det)
            tr <- xml_add_child(x, "trans", from = "D", to = "R", type = "gamma", mean = "recovery_period", shape = shape$rec)
            if (time_step == 0) {
                pass_events <- clamp(pass_events, 2, 3)
                if (pass_events >= 3) xml_set_attr(ts, "data_column", pass_events + 1)
                xml_set_attr(td, "data_column", pass_events + 2)
                xml_set_attr(tr, "data_column", pass_events + 3)
            }
        },
        "SEIDR" = {
            ts <- xml_add_child(x, "trans", from = "S", to = "E", type = "infection", beta = "beta", inf_model = "frequency dependent")
            tl <- xml_add_child(x, "trans", from = "E", to = "I", type = "gamma", mean = "latent_period", shape = shape$lat)
            td <- xml_add_child(x, "trans", from = "I", to = "D", type = "gamma", mean = "detection_period", shape = shape$det)
            tr <- xml_add_child(x, "trans", from = "D", to = "R", type = "gamma", mean = "recovery_period", shape = shape$rec)
            if (time_step == 0) {
                pass_events <- clamp(pass_events, 2, 4)
                if (pass_events >= 4) xml_set_attr(ts, "data_column", pass_events)
                if (pass_events >= 3) xml_set_attr(tl, "data_column", pass_events + 1)
                xml_set_attr(td, "data_column", pass_events + 2)
                xml_set_attr(tr, "data_column", pass_events + 3)
            }
        }
    )
    
    ### Add individual effects ----
    slidr <- c("s", "l", "i", "d", "r")
    slidr_ie <- paste0(slidr, "_a,", slidr, "_e")
    names(slidr_ie) <- slidr
    
    if ("susceptibility" %in% traitnames) {
        xml_set_attr(ts, "individual_effect", slidr_ie[["s"]])
    }
    
    if ("infectivity" %in% traitnames) {
        xml_set_attr(ci, "individual_effect", slidr_ie[["i"]])
        if (exists("cd")) {
            xml_set_attr(cd, "individual_effect", slidr_ie[["i"]])
        }
    }
    
    if ("recoverability" %in% traitnames) {
        xml_set_attr(tr, "individual_effect", slidr_ie[["r"]])
        
        # link latency and detectability
        if ("latency" %in% traitnames && substr(link_traits, 2, 2) == "r") {
            xml_set_attr(tl, "individual_effect", slidr_ie[["r"]])
        }
        
        if ("detectability" %in% traitnames && substr(link_traits, 4, 4) == "r") {
            xml_set_attr(td, "individual_effect", slidr_ie[["r"]])
        }
    }
    
    
    
    
    ### Add fixed effects ----
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
    if (fex != "") {
        xml_set_attr(ts, "fixed_effect", fex)
    }
    
    fex <- fe_str("l")
    if (fex != "") {
        xml_set_attr(tl, "fixed_effect", fex)
    }
    
    fex <- fe_str("i")
    if (fex != "") {
        xml_set_attr(ci, "fixed_effect", fex)
        if (exists("cd")) {
            xml_set_attr(cd, "fixed_effect", fex)
        }
    }
    
    fex <- fe_str("d")
    if (fex != "") {
        xml_set_attr(td, "fixed_effect", fex)
    }
    
    fex <- fe_str("r")
    if (fex != "") {
        xml_set_attr(tr, "fixed_effect", fex)
    }
    
    
    ### Group effect ----
    if (group_effect < 0) {
        # useful if group_effect has been set outside of make_parameters()
        priors[parameter == "sigma", use := FALSE]
    } else {
        xml_add_child(x, xml_comment("Group effect"))
        xml_add_child(x, "group_effect", sigma = "sigma")
    }
    
    
    ## Covariance ----
    
    # Note that we have two child nodes both called "covariance", so
    # xml_add_child() fails when we add to the second one. Instead we use
    # the position which is also the length of `x`.
    
    
    sir <- c("susceptibility", "infectivity", "recoverability")
    used_traits <- sir %in% traitnames
    n_used <- sum(used_traits)
    use_cov <- any(used_traits)
    
    if (use_cov) {
        # generate a string like "s_a,i_a" depending on which traits are used
        ie_a_str <- paste0(c("s", "i", "r")[used_traits], "_a", collapse = ",")
        ie_e_str <- gsub("a", "e", ie_a_str)
        
        cov_G_mat <- matrix(paste0("cov_G_", c("ss", "ii", "rr")[used_traits]), ncol = 1)
        cov_E_mat <- sub("G", "E", cov_G_mat)
        
        cor_G_mat <- diag(1, n_used, n_used)
        cors <- c(if (all(used_traits[c(1, 2)])) "r_G_si",
                  if (all(used_traits[c(1, 3)])) "r_G_sr",
                  if (all(used_traits[c(2, 3)])) "r_G_ir")
        cor_G_mat[upper.tri(cor_G_mat)] <- cors
        cor_G_mat[lower.tri(cor_G_mat)] <- t(cor_G_mat)[lower.tri(cor_G_mat)]
        cor_E_mat <- sub("G", "E", cor_G_mat)
        
        
        # Genetic Covariance
        xml_add_child(x, xml_comment("Genetic covariance between different individual effects"))
        xml_add_child(x, "covariance",
                      individual_effect = ie_a_str,
                      relationship_matrix = "A")
        xml_add_child(xml_child(x, length(xml_children(x))),
                      "variance",
                      table_to_tsv_string(cov_G_mat))
        xml_add_child(xml_child(x, length(xml_children(x))),
                      "correlation",
                      table_to_tsv_string(cor_G_mat))
        
        # Environmental Covariance
        xml_add_child(x, xml_comment("Environmental covariance between different individual effects"))
        xml_add_child(x, "covariance",
                      individual_effect = ie_e_str,
                      relationship_matrix = "I")
        xml_add_child(xml_child(x, length(xml_children(x))),
                      "variance",
                      table_to_tsv_string(cov_E_mat))
        xml_add_child(xml_child(x, length(xml_children(x))),
                      "correlation",
                      table_to_tsv_string(cor_E_mat))
    }
    
    ## Demography ----
    xml_add_child(x, xml_comment("No. of sires + progeny, no. of groups"))
    xml_add_child(x, "data", N = nsires + nprogeny)
    xml_add_child(x, "data", Z = ngroups)
    
    
    ## Observation and Inference periods ----
    xml_add_child(x, xml_comment("Inference and observation periods"))
    if (use_fb_data) {
        if (setup == "fishboost") {
            for (i in 1:72) {
                tmax_val <- if (i <= 36) 104 else 160
                xml_add_child(x, "inference", group = i, tmin = 0, tmax = tmax_val)
                xml_add_child(x, "observation", group = i, tmin = 0, tmax = tmax_val)
            }
        } else if (setup == "fb1") {
            xml_add_child(x, "inference", tmin = 0, tmax = 104)
            xml_add_child(x, "observation", tmin = 0, tmax = 104)
        } else {
            xml_add_child(x, "inference", tmin = 0, tmax = 160)
            xml_add_child(x, "observation", tmin = 0, tmax = 160)
        }
    } else {
        xml_add_child(x, "inference", tmin = 0, tmax = "infinity")
        xml_add_child(x, "observation", tmin = 0, tmax = "infinity")
    }
    
    
    
    ## Priors ----
    xml_add_child(x, xml_comment("Model priors"))
    for (i in seq_len(nrow(priors))) {
        ppi <- as.list(priors[i, ])
        
        # skip over unused parameters
        if (ppi$use == FALSE) next
        
        xml_add_child(x, "prior",
                      parameter = ppi$parameter,
                      type = "Flat", val1 = ppi$val1, val2 = ppi$val2)
    }
    
    
    ## Event times ----
    
    xml_add_child(x, xml_comment("Information about individuals"))

    if (time_step == 0) {
        xdt <- xml_add_child(x, "datatable",
                             id = 1, group = 2, initial_comp = 3, #s_a = 6, i_a = 7, r_a = 8,
                             "\nREPLACE_WITH_DATA\n")
    } else {
        xdt <- xml_add_child(x, "datatable",
                             id = 1, group = 2, initial_comp = 3, comp_status = 4,
                             #s_a = 5, i_a = 6, r_a = 7,
                             "\nREPLACE_WITH_DATA\n")
    }
    
    
    ## DT col names ----
    if (!use_fb_data) {
        if (used_traits[1]) xml_set_attr(xdt, "s_a", which(names(data) == "susceptibility_BV"))
        if (used_traits[2]) xml_set_attr(xdt, "i_a", which(names(data) == "infectivity_BV"))
        if (used_traits[3]) xml_set_attr(xdt, "r_a", which(names(data) == "recoverability_BV"))
    }
    
    # add <datatable ... trial_r="9", donor_i="10"> etc.
    
    
    for (str in unique(trial_names[strsplit(trial_fe, "")[[1]]])) {
        xml_set_attr(xdt, str, which(names(data) == "trial_fe"))
    }
    
    for (str in unique(donor_names[strsplit(donor_fe, "")[[1]]])) {
        xml_set_attr(xdt, str, which(names(data) == "donor_fe"))
    }
    
    for (str in unique(txd_names[strsplit(txd_fe, "")[[1]]])) {
        xml_set_attr(xdt, str, which(names(data) == "txd_fe"))
    }
    
    # for (i in strsplit(trial_fe, "")[[1]])
    # xml_set_attr(xdt, paste0("trial_", i), which(names(data) == "trial_fe"))
    # for (i in strsplit(donor_fe, "")[[1]])
    # xml_set_attr(xdt, paste0("donor_", i), which(names(data) == "donor_fe"))
    # for (i in strsplit(txd_fe, "")[[1]])
    # xml_set_attr(xdt, paste0("txd_", i), which(names(data) == "txd_fe"))
    
    
    # Write the data to a TSV file, then use sed to replace the string
    # "REPLACE_WITH_DATA". This is trickier, but *much* faster than trying
    # to generate a huge string in the XML stream.
    
    data_file <- glue("{data_dir}/{name}-data.tsv")
    write.table(data, file = data_file, sep = "\t", quote = FALSE,
                row.names = FALSE, col.names = FALSE)
    
    
    ## GRM ----
    
    # Clip out the dams
    nondams <- c(seq.int(nsires), seq.int(nparents + 1L, ntotal))
    GRM_nd <- GRM[nondams, nondams]
    
    # Convert the GRM into a sparse matrix, use `summary` to extract values,
    # then convert that into a data table
    
    GRM_nd_dt <- as.data.table(summary(Matrix(GRM_nd, sparse = TRUE)))
    
    # Use the DT to make 0-based indices
    GRM_nd_dt[, `:=`(i = i - 1L, j = j - 1L)]
    
    # Account for summary only giving UT matrix, and SIRE wanting access to
    # both the LT and UT parts
    GRM_nd_dt <- rbind(GRM_nd_dt, GRM_nd_dt[i != j, .(i = j, j = i, x)])
    setkey(GRM_nd_dt, j, i)
    
    
    # write the DT as a TSV file
    A_file <- glue("{data_dir}/{name}-A.tsv")
    write.table(GRM_nd_dt, file = A_file, sep = "\t", quote = FALSE,
                row.names = FALSE, col.names = FALSE)
    
    xml_add_child(x, xml_comment("A-matrix or GRM"))
    xml_add_child(x, "A_nonzero", name = "A", "\nREPLACE_WITH_A\n")
    
    
    ## Prediction Accuracies ----
    
    
    # build up ind = "Ind0,Ind1,...,Ind99"
    sire_inds <- paste0("Ind", seq.int(0L, nsires - 1L), collapse = ",")
    xml_add_child(x, xml_comment("Prediction Accuracies for sires"))
    xml_add_child(x, "prediction_accuracy", name = "sire", ind = sire_inds)
    
    # build up ind = "Ind100,Ind101,...,Ind2099"
    prog_inds <- paste0("Ind", seq.int(nsires, nsires + nprogeny - 1L), collapse = ",")
    xml_add_child(x, xml_comment("Prediction Accuracies for progeny"))
    xml_add_child(x, "prediction_accuracy", name = "progeny", ind = prog_inds)
    
    
    # Write the XML file ----
    xml_file <- glue("{data_dir}/{name}.xml")
    write_xml(x, xml_file) #options = c("as_xml", "format") #?
    
    
    ## Inject the TSV files using sed ----
    cmd1 <- paste0("sed -e '/REPLACE_WITH_DATA/ {' -e 'r ", data_file, "' -e  'd' -e '}' -i ", xml_file)
    cmd2 <- paste0("sed -e '/REPLACE_WITH_A/ {' -e 'r ", A_file, "' -e  'd' -e '}' -i ", xml_file)
    
    # Remove indents (they're terrible anyway), and insert newlines before
    # comments for readability
    cmd3 <- paste0("sed -i 's/^\\t*//; s/^ *//; s/\\(<!--\\)/\\n\\1/' ", xml_file)
    
    # On Mac need to use GNU sed (gsed) instead of BSD sed
    if (Sys.info()["sysname"] == "Darwin") {
        cmd1 <- glue("g{cmd1}")
        cmd2 <- glue("g{cmd2}")
        cmd3 <- glue("g{cmd3}")
    }
    
    system(cmd1)
    system(cmd2)
    system(cmd3)
    
    ## Remove temporary files ----
    system(glue("rm {A_file} {data_file}"))
    
    x
}
