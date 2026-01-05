{
    library(data.table)
    library(stringr)
}

# Potentially reclassify donors as recipients
fix_fb_data <- function(popn, params) {
    {
        fix_donors   <- params$fix_donors
        fix_eq_time  <- params$fix_eq_time
        t_demote     <- params$t_demote
        tmax         <- params$tmax
        compartments <- params$compartments
    }

    message("Correcting FB data ...")
    
    # Check in case passing old params
    t_demote <- t_demote %||% c(Inf, Inf)
    fix_eq_time <- fix_eq_time %||% FALSE
    if ("use_parasites" %in% params && params$use_parasites != "") fix_donors <- "no_Tsym_survivors"
    
    popn2 <- copy(popn)
    
    # Remove deaths at tmax
    message(" - Fixing Tdeath >= Tmax ...")
    popn2[Tdeath >= tmax[str_glue("t{trial}")], Tdeath := NA]
    
    # Fix Tsym = Tdeath by putting Tsym down 1/2 a day
    if (fix_eq_time) {
        message(" - Fixing Tsym = Tdeath ...")
        popn2[Tsym == Tdeath, Tsym := Tsym - 0.5]
    }
    
    
    if ("time" %in% fix_donors) {
        message(" - Demoting donors with Tsym > ", str_flatten_comma(t_demote))
        popn2[donor == 1 & Tsym > t_demote[trial], `:=`(donor = 0, Tinf = NA)]
    }
    
    new_val <- if ("set_to_R" %in% fix_donors) length(compartments) -1 else 0
    
    if ("no_Tsym_survivors" %in% fix_donors) {
        message(" - Demoting surviving donors with no symptoms ...")
        popn2[donor == 1 & is.na(Tsym) & is.na(Tdeath), `:=`(donor = new_val, Tinf = NA)]
    }
    
    # Keep at least 1 donor from each group, even if none show parasites
    zero_groups <- popn2[, .(Nd = sum(donor == 1)), group][Nd == 0, group]
    if (length(zero_groups) > 0) {
        message(" - Groups with no donors: ", str_flatten_comma(zero_groups))
        ids_to_keep <- popn2[group %in% zero_groups, id[1], group][, V1]
        popn2[ids_to_keep, `:=`(donor = 1, Tinf = 0)]
    }
    
    if ("set_to_R" %in% fix_donors) {
        popn2[donor == new_val, Tdeath := NA]
    }
    
    popn2
}
