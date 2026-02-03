{
    library(stringr)
}

remove_bici_fes <- function(params) {
    trial_fe <- params$trial_fe
    donor_fe <- params$donor_fe
    
    params2 <- copy(params)
    priors <- params2$priors
    
    message("Checking if parameters need merged for BICI if no Trial or Donor FEs")
    
    # Trial ----
    if (params$ntrials > 1) {
        if (!str_detect(trial_fe, "s|i")) {
            message("- Merging beta across trials")
            priors[str_detect(parameter, "beta_Tr"), true_val := mean(true_val)]
        }
        if (!str_detect(trial_fe, "l")) {
            message("- Merging LP across trials")
            priors[str_detect(parameter, "lat.*Don"), true_val := mean(true_val)]
            priors[str_detect(parameter, "lat.*Rec"), true_val := mean(true_val)]
        }
        if (!str_detect(trial_fe, "d")) {
            message("- Merging DP across trials")
            priors[str_detect(parameter, "det.*Don"), true_val := mean(true_val)]
            priors[str_detect(parameter, "det.*Rec"), true_val := mean(true_val)]
        }
        if (!str_detect(trial_fe, "t")) {
            message("- Merging RP across trials")
            priors[str_detect(parameter, "rem.*Don"), true_val := mean(true_val)]
            priors[str_detect(parameter, "rem.*Rec"), true_val := mean(true_val)]
        }
    }
        
    # Donor ----
    
    # Need to weight by proportion of inoculated status)
    pD <- params$I0 / params$group_size
    wtI <- c(pD, 1 - pD)
    
    if (!str_detect(donor_fe, "l")) {
        message("- Merging LP across inoculation status")
        priors[str_detect(parameter, "lat.*Tr1"), true_val := sum(wtI * true_val)]
        priors[str_detect(parameter, "lat.*Tr2"), true_val := sum(wtI * true_val)]
    }
    if (!str_detect(donor_fe, "d")) {
        message("- Merging DP across inoculation status")
        priors[str_detect(parameter, "det.*Tr1"), true_val := sum(wtI * true_val)]
        priors[str_detect(parameter, "det.*Tr2"), true_val := sum(wtI * true_val)]
    }
    if (!str_detect(donor_fe, "t")) {
        message("- Merging RP across inoculation status")
        priors[str_detect(parameter, "rem.*Tr1"), true_val := sum(wtI * true_val)]
        priors[str_detect(parameter, "rem.*Tr2"), true_val := sum(wtI * true_val)]
    }
    
    params2
}
