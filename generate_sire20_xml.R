# if we want the same data as example.xml (handy for exploring the structure)

# sire_ex <- read_xml("example.xml")
# sire_ex_list <- as_list(sire_ex)

generate_sire20_xml <- function(data, GRM, params) {
    name       <- params$name
    nsires     <- params$nsires
    nprogeny   <- params$nprogeny
    nparents   <- params$nparents
    ntotal     <- params$ntotal
    ngroups    <- params$ngroups
    model_type <- params$model_type
    nsample    <- params$nsample
    burnin     <- params$burnin
    t_gap      <- params$t_gap
    data_dir   <- params$data_dir

    message(glue("Generating SIRE 2.0 XML file {data_dir}/{name}.xml ..."))

    # Build XML file `x` ----
    x <- xml_new_root("SIRE", version = "2.0")


    ## MCMC options ----
    xml_add_child(x, xml_comment("MCMC options"))
    xml_add_child(x, "mcmc", nsample = nsample, burnin = burnin)


    ## Model details ----

    # SIRE 2.0 only supports SIS and SIR models
    if (!(model_type %in% c("SIR", "SIS"))) {
        message(" - pretending this is an SIR model")
        model_type <- "SIR"
    }

    xml_add_child(x, xml_comment("Model details"))
    xml_add_child(x, "model", type = model_type, residual = "on", groupeff = "on")
    xml_add_child(x, "data", N = nsires + nprogeny)
    xml_add_child(x, "data", Z = ngroups)

    # recover tmax from strings
    tmax <- data[!(Trec %in% c(".",  "no", "NA")),  ceiling(max(as.numeric(Trec)))]

    xml_add_child(x, "inference", tmin = 0, tmax = tmax)
    xml_add_child(x, "observation", tmin = 0, tmax = tmax)


    ## Priors ----

    # SIRE 2.0 expects gamma_shape to be called "k"
    params$priors[parameter == "gamma_shape", symbol := "k"]

    xml_add_child(x, xml_comment("Model priors"))
    for (i in seq_len(nrow(params$priors))) {
        ppi <- as.list(params$priors[i, ])

        # skip over parameters SIRE 2.0 doesn't recognise
        if (ppi$use_s20 == FALSE) next

        xml_add_child(x, "prior",
                      parameter = ppi$symbol,
                      type = "Flat", val1 = ppi$val1, val2 = ppi$val2)
    }


    ## Event times ----

    xml_add_child(x, xml_comment("Information about individuals"))
    # if t_gap > 0 then It/Rt need to be state
    if (t_gap > 0) {
        xml_add_child(x, "datatable",
                      id = 1, group = 2, state = 3, qg = 4, qf = 5, qr = 6,
                      "\nREPLACE_WITH_DATA\n")
    } else {
        xml_add_child(x, "datatable",
                      id = 1, group = 2, It = 3, Rt = 4, qg = 5, qf = 6, qr = 7,
                      "\nREPLACE_WITH_DATA\n")
    }

    data_file <- glue("{data_dir}/{name}-data.tsv")
    write.table(data, file = data_file, sep = "\t", quote = FALSE,
                row.names = FALSE, col.names = FALSE)

    ## Relationship Matrix ----

    # clip out the dams
    nondams <- c(seq.int(nsires), seq.int(nparents + 1L, ntotal))
    GRM_nd <- GRM[nondams, nondams]

    # convert the GRM into a sparse matrix, use `summary` to extract values,
    # then convert that into a data table
    GRM_nd_dt <- as.data.table(summary(Matrix(GRM_nd, sparse = TRUE)))

    # use the DT to fix the indices
    GRM_nd_dt[, `:=`(i = i - 1L, j = j - 1L)]

    # write the DT as a TSV file
    A_file <- glue("{data_dir}/{name}-A.tsv")
    write.table(GRM_nd_dt, file = A_file, sep = "\t", quote = FALSE,
                row.names = FALSE, col.names = FALSE)

    xml_add_child(x, xml_comment("A-matrix or GRM"))
    xml_add_child(x, "A_nonzero", "\nREPLACE_WITH_A\n")


    ## Prediction Accuracies ----

    # build up ind = "Ind0,Ind1,...,Ind99"
    sire_inds <- paste0("Ind", seq.int(0L, nsires - 1L), collapse = ",")
    xml_add_child(x, xml_comment("Prediction Accuracies for sires"))
    xml_add_child(x, "prediction_accuracy", name = "sire", ind = sire_inds)

    # build up ind = "Ind100,Ind101,...,Ind2099"
    prog_inds <- paste0("Ind", seq.int(nsires, nsires + nprogeny - 1L), collapse = ",")
    xml_add_child(x, xml_comment("Prediction Accuracies for progeny"))
    xml_add_child(x, "prediction_accuracy", name = "progeny", ind = prog_inds)


    # Write XML file ----
    xml_file <- glue("{data_dir}/{name}.xml")
    write_xml(x, xml_file)


    ## Inject the TSV files with sed ----
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

    # Tidy up the remaining TSV files
    system(glue("rm {A_file} {data_file}"))

    x
}
