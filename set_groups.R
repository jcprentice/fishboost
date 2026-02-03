# Copy the popn DT and assign groups, donors and recipients, group effect, and
# trials.

set_groups <- function(popn, params) {
    message("Assigning groups to progeny ...")
    
    {
        setup        <- params$setup
        nprogeny     <- params$nprogeny
        ngroups      <- params$ngroups
        group_layout <- params$group_layout
        group_size   <- params$group_size
        group_effect <- params$group_effect
        I0           <- params$I0
        msgs         <- params$msgs
    }
    
    # Group layout will break if set to fishboost unless the pedigree matches
    if (group_layout == "fishboost" && !str_starts(setup, "fb")) {
        group_layout <- "random"
    }
    
    # New DT with just progeny IDs
    dt <- popn[sdp == "progeny", .(id)]
    
    # Set trial
    if (str_starts(setup, "fb")) {
        fb_data <- readRDS(str_glue("fb_data/{setup}.rds"))
        dt[, trial := fb_data[dt$id, trial]]
    } else {
        dt[, trial := 1L]
    }
    
    # Set groups
    switch(group_layout,
           "random" = {
               if (msgs) message("- shuffling groups")
               dt[, group := sample(rep(1:ngroups, length.out = .N))]
           }, "family" = {
               if (msgs) message("- by family")
               dt[, group := rep(1:ngroups, each = group_size, length.out = .N)]
           }, "striped" = {
               if (msgs) message("- by stripes")
               dt[, group := rep(1:ngroups, length.out = .N)]
           }, "fishboost" = {
               if (msgs) message("- copying Fishboost groups")
               dt[, group := fb_data[dt$id, group]]
           }, stop("Unrecognised group_layout!")
    )
    
    # Assign I0 initial infectives to each group
    dt[, donor := {
        x = integer(.N)
        x[1:min(I0, .N)] = 1L
        x
    }, by = group]
    
    # Handle group effect
    # note: group_effect is a multiplier, so x1 means no effect
    if (group_effect >= 0) {
        # also need to handle the case where groups are not of equal size
        if (msgs) message("- group_effect = ", group_effect)
        dt[, GE := exp(rnorm(1L, 0, group_effect)), by = group]
    } else {
        dt[, GE := 1]
    }
    
    # Copy popn to popn2 by merging on "id" and excluding columns in common
    # (which will be overwritten from dt)
    cols <- c("id", setdiff(names(popn), names(dt)))
    popn2 <- merge(popn[, ..cols], dt, by = "id", all = TRUE)
    
    popn2
}


