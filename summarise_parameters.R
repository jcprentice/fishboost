
summarise_parameters <- function(params) {
    with(params, {
        if (sim_new_data == "no") {
            message("Using Fishboost data")
        } else if (sim_new_data == "summary_inf") {
            message("Using data simulated in BICI from prior")
        } else if (sim_new_data == "summary_ps") {
            message("Using data simulated in BICI from posterior")
        } else {
            message(str_glue("Simulating new {model_type} data via {x}",
                             x = str_to_upper(sim_new_data)))
        }

        message(str_glue(
            "- Demography is:",
            "\t{nsires} sires, {ndams} dams, {nprogeny} progeny ({ntotal} total)",
            "\t{ngroups} group{s} (group size {group_size})",
            s = if (ngroups > 1) "s" else "",
            .trim = FALSE, .sep = "\n"
        ))

        message(str_glue(
            "- Individual Effects are:",
            "\t'{ts}' -> '{ies}'",
            ts = str_flatten(names(model_traits)),
            ies = if (IEs %in% c("", "none")) "none" else IEs,
            .trim = FALSE, .sep = "\n"
        ))

        message("- Sigma_G: ")
        capture_message(Sigma_G[model_traits, model_traits])

        message(str_glue(
            "- Fixed Effects are:",
            "\tTrial = '{tfe}', Donor = '{dfe}', TxD = '{xfe}', Weight = '{wfe}'",
            tfe = if (FEs$trial  %in% c("", "none")) "none" else FEs$trial,
            dfe = if (FEs$donor  %in% c("", "none")) "none" else FEs$donor,
            xfe = if (FEs$txd    %in% c("", "none")) "none" else FEs$txd,
            wfe = if (FEs$weight %in% c("", "none")) "none" else FEs$weight,
            .trim = FALSE, .sep = "\n"
        ))

        message(str_glue(
            "- Estimated R0:",
            "\t{r0}",
            r0 = signif(R0, 3),
            .trim = FALSE, .sep = "\n"
        ))

        if (patch_dataset != "") {
            message(str_glue(
                "- Patching with dataset {patch_dataset} / {patch_name}",
                "\ttype = {patch_type}, state = {patch_state}",
                .trim = FALSE, .sep = "\n"
            ))
        }

        if (bici_cmd == "inf") {
            message(str_glue(
                "- Running MCMC with:",
                "\t{ns} updates, {th} samples, {burnprop} burnin, {nchains} chains",
                "- BICI script file:",
                "\t'{data_dir}/{name}.bici'",
                "- Results file:",
                "\t'{results_dir}/{name}.rds'",
                ns = format(nsample, scientific = FALSE, big.mark = ","),
                th = format(thinto,  scientific = FALSE, big.mark = ","),
                .trim = FALSE, .sep = "\n"
            ))
        }
    })
}

