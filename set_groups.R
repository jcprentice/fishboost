# Copy the traits DT and assign groups, donors and recipients, group effect, and
# trials.

set_groups <- function(traits, params) {
    message("Assigning groups to progeny")

    {
        setup        <- params$setup
        nprogeny     <- params$nprogeny
        ngroups      <- params$ngroups
        group_layout <- params$group_layout
        group_size   <- params$group_size
        group_effect <- params$group_effect
        I0           <- params$I0
    }

    # If using a Fishboost layout...
    using_fb_setup <- setup %in% c("fishboost", "fb1", "fb2")
    if (using_fb_setup) {
        fb_data <- switch(
            setup,
            "fishboost" = readRDS("fb_data/fb_data12.rds"),
            "fb1"       = readRDS("fb_data/fb_data1.rds"),
            "fb2"       = readRDS("fb_data/fb_data2.rds")
        )
    }


    # Group layout will break if set to fishboost unless the pedigree matches
    if (group_layout == "fishboost" && !using_fb_setup) {
        group_layout <- "random"
    }

    # New DT with just progeny IDs
    dt <- traits[sdp == "progeny", .(id)]

    # Set trial
    if (using_fb_setup) {
        dt[, trial := fb_data[dt$id, trial]]
    } else {
        dt[, trial := 1L]
    }

    # Set groups
    switch(group_layout,
           "random" = {
               message(" - shuffling groups")
               dt[, group := sample(rep(1:ngroups, length.out = .N))]
           }, "family" = {
               message(" - by family")
               dt[, group := rep(1:ngroups, each = group_size, length.out = .N)]
           }, "striped" = {
               message(" - by stripes")
               dt[, group := rep(1:ngroups, length.out = .N)]
           }, "fishboost" = {
               message(" - copying Fishboost groups")
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
        message(" - group_effect = ", group_effect)
        dt[, GE := exp(rnorm(1L, 0, group_effect)), by = group]
    } else {
        dt[, GE := 1]
    }

    # Copy traits to pop by merging on "id" and excluding columns in common
    # (which will be overwritten from dt)
    cols <- c("id", setdiff(names(traits), names(dt)))
    pop <- merge(traits[, ..cols], dt, by = "id", all = TRUE)

    pop
}


