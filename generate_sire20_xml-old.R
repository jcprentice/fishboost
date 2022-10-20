# if we want the same data as example.xml (handy for exploring the structure)

# sire_ex <- read_xml("example.xml")
# sire_ex_list <- as_list(sire_ex)

generate_sire20_xml <- function(data, GRM, params) {
    # name          <- params$name
    # nsires        <- params$nsires
    # nprogeny      <- params$nprogeny
    # nparents      <- params$nparents
    # ngroups       <- params$ngroups
    # model_type    <- params$model_type
    # nsample       <- params$nsample
    # burnin        <- params$burnin
    # t_gap         <- params$t_gap
    
    with (params, {
        message("Generating XML file data/", name, ".xml ...")
        
        # recover tmax from strings
        tmax <- data[!(Trec %in% c(".", "no")), ceiling(max(as.numeric(Trec)))]
        
        # Turn DT into an XML node ----
        
        if (t_gap > 0) {
            data_xml <- structure(
                list(data_table_to_tsv_string(data)),
                id = 1, group = 2, state = 3, qg = 4, qf = 5, qr = 6)
        } else {
            data_xml <- structure(
                list(data_table_to_tsv_string(data)),
                id = 1, group = 2, It = 3, Rt = 4, qg = 5, qf = 6, qr = 7)
        }
        
        # GRM ----
        
        # clip out the dams
        nondams <- c(seq.int(nsires), seq.int(nprogeny) + nparents)
        GRM2 <- GRM[nondams, nondams]
        
        # turn the GRM into a sparse matrix
        sA <- summary(Matrix(GRM2, sparse = TRUE))
        
        # then turn that into a matrix that we can write
        mA <- matrix(c(sA$i - 1, sA$j - 1, sA$x), ncol = 3)
        
        # turn this into an XML node
        A_nonzero_xml <- structure(list(data_table_to_tsv_string(mA)))
        
        
        # Prediction Accuracies ----
        
        # build up ind = "Ind0,Ind1,...,Ind99"
        sire_vals <- seq.int(0L, nsires - 1L)
        sire_inds <- paste0("Ind", sire_vals, collapse = ",")
        pa_sire_xml <- structure(list(), name = "sire", ind = sire_inds)
        
        # build up ind = "Ind100,Ind101,...,Ind2099"
        prog_vals <- seq.int(nsires, nsires + nprogeny - 1L)
        prog_inds <- paste0("Ind", prog_vals, collapse = ",")
        pa_progeny_xml <- structure(list(), name = "progeny", ind = prog_inds)
        
        
        # SIRE 2.0 only supports SIS and SIR models
        if (!(model_type %in% c("SIR", "SIS"))) {
            message(" - pretending this is an SIR model")
            model_type <- "SIR"
        }
        
        x <- list(
            SIRE = structure(
                version = "2.0",
                list(
                    mcmc = structure(list(), nsample = nsample, burnin = burnin),
                    model = structure(list(), type = model_type, residual = "on", groupeff = "on"),
                    data = structure(list(), N = nsires + nprogeny),
                    data = structure(list(), Z = ngroups),
                    inference = structure(list(), tmin = 0, tmax = tmax),
                    observation = structure(list(), tmin = 0, tmax = tmax),
                    
                    prior = structure(list(), parameter = "β",    type = "Flat", val1 = prior_beta1,  val2 = prior_beta2),
                    prior = structure(list(), parameter = "γ",    type = "Flat", val1 = prior_gamma1, val2 = prior_gamma2),
                    prior = structure(list(), parameter = "k",    type = "Flat", val1 = prior_k1,     val2 = prior_k2),
                    prior = structure(list(), parameter = "G",    type = "Flat", val1 = -3, val2 = 3),
                    prior = structure(list(), parameter = "σ_G",  type = "Flat", val1 = 0,  val2 = 0.5),
                    prior = structure(list(), parameter = "q_g",  type = "Flat", val1 = -5, val2 = 5),
                    prior = structure(list(), parameter = "q_f",  type = "Flat", val1 = -5, val2 = 5),
                    prior = structure(list(), parameter = "q_r",  type = "Flat", val1 = -5, val2 = 5),
                    prior = structure(list(), parameter = "Ω_gg", type = "Flat", val1 = prior_Sigma_gg1, val2 = prior_Sigma_gg2),
                    prior = structure(list(), parameter = "Ω_ff", type = "Flat", val1 = prior_Sigma_ff1, val2 = prior_Sigma_ff2),
                    prior = structure(list(), parameter = "Ω_rr", type = "Flat", val1 = prior_Sigma_rr1, val2 = prior_Sigma_rr2),
                    prior = structure(list(), parameter = "Ω_gf", type = "Flat", val1 = prior_Sigma_gf1, val2 = prior_Sigma_gf2),
                    prior = structure(list(), parameter = "Ω_gr", type = "Flat", val1 = prior_Sigma_gr1, val2 = prior_Sigma_gr2),
                    prior = structure(list(), parameter = "Ω_fr", type = "Flat", val1 = prior_Sigma_fr1, val2 = prior_Sigma_fr2),
                    prior = structure(list(), parameter = "ε_g",  type = "Flat", val1 = -5, val2 = 5),
                    prior = structure(list(), parameter = "ε_f",  type = "Flat", val1 = -5, val2 = 5),
                    prior = structure(list(), parameter = "ε_r",  type = "Flat", val1 = -5, val2 = 5),
                    prior = structure(list(), parameter = "Ψ_gg", type = "Flat", val1 = prior_Psi_gg1, val2 = prior_Psi_gg2),
                    prior = structure(list(), parameter = "Ψ_ff", type = "Flat", val1 = prior_Psi_ff1, val2 = prior_Psi_ff2),
                    prior = structure(list(), parameter = "Ψ_rr", type = "Flat", val1 = prior_Psi_rr1, val2 = prior_Psi_rr2),
                    prior = structure(list(), parameter = "Ψ_gf", type = "Flat", val1 = prior_Psi_gf1, val2 = prior_Psi_gf2),
                    prior = structure(list(), parameter = "Ψ_gr", type = "Flat", val1 = prior_Psi_gr1, val2 = prior_Psi_gr2),
                    prior = structure(list(), parameter = "Ψ_fr", type = "Flat", val1 = prior_Psi_fr1, val2 = prior_Psi_fr2),
                    
                    datatable = data_xml,
                    # datatable = dt_chris_xml, # from Chris's example.xml
                    # Ainv_nonzero = Ainv_chris_xml,
                    A_nonzero = A_nonzero_xml,
                    prediction_accuracy = pa_sire_xml,
                    prediction_accuracy = pa_progeny_xml
                )
            )
        )
    })
    
    
    
    xml_x <- as_xml_document(x)
    write_xml(xml_x, paste0("data/", name, ".xml"), options = c("as_xml", "format"))
    
    xml_x
}
