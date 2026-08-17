{
    library(stringr)
}

remove_bici_fes <- function(params) {
    FEs <- params$FEs
    trial <- if ("trial" %in% names(FEs)) FEs$trial else ""
    donor <- if ("donor" %in% names(FEs)) FEs$donor else ""

    params2 <- copy(params)
    priors <- params2$priors

    message("Checking if parameters need merged for BICI if no Trial or Donor FEs")

    # Trial ----
    if (params$ntrials > 1) {
        if (!str_detect(trial, "[si]")) {
            message("- Merging beta across trials")
            priors[str_detect(parameter, "beta_Tr"), true_val := mean(true_val)]
        }
        if (!str_detect(trial, "l")) {
            message("- Merging LP across trials")
            priors[str_detect(parameter, "lat.*Don"), true_val := mean(true_val)]
            priors[str_detect(parameter, "lat.*Rec"), true_val := mean(true_val)]
        }
        if (!str_detect(trial, "d")) {
            message("- Merging DP across trials")
            priors[str_detect(parameter, "det.*Don"), true_val := mean(true_val)]
            priors[str_detect(parameter, "det.*Rec"), true_val := mean(true_val)]
        }
        if (!str_detect(trial, "t")) {
            message("- Merging RP across trials")
            priors[str_detect(parameter, "rem.*Don"), true_val := mean(true_val)]
            priors[str_detect(parameter, "rem.*Rec"), true_val := mean(true_val)]
        }
    }

    # Donor ----

    # Need to weight by proportion of inoculated status)
    pD <- params$I0 / params$group_size
    wtI <- c(pD, 1 - pD)

    if (!str_detect(donor, "l")) {
        message("- Merging LP across inoculation status")
        priors[str_detect(parameter, "LP.*Tr1"), true_val := sum(wtI * true_val)]
        priors[str_detect(parameter, "LP.*Tr2"), true_val := sum(wtI * true_val)]
    }
    if (!str_detect(donor, "d")) {
        message("- Merging DP across inoculation status")
        priors[str_detect(parameter, "DP.*Tr1"), true_val := sum(wtI * true_val)]
        priors[str_detect(parameter, "DP.*Tr2"), true_val := sum(wtI * true_val)]
    }
    if (!str_detect(donor, "t")) {
        message("- Merging RP across inoculation status")
        priors[str_detect(parameter, "RP.*Tr1"), true_val := sum(wtI * true_val)]
        priors[str_detect(parameter, "RP.*Tr2"), true_val := sum(wtI * true_val)]
    }

    params2
}

