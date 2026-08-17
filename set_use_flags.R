safe_source("utils.R")

set_use_flags <- function(params) {
    {
        all_traits    <- params$all_traits
        model_traits  <- params$model_traits
        IEs           <- params$IEs
        FEs           <- params$FEs
        setup         <- params$setup
        LP_dist       <- params$LP_dist
        DP_dist       <- params$DP_dist
        RP_dist       <- params$RP_dist
        group_effect  <- params$group_effect
        # DTs are modified by reference, so this doesn't need to be copied back
        priors        <- params$priors
    }

    message("Setting use flags in priors ...")

    priors[, use := FALSE]

    # Turn on base parameters

    # get letters for each trait
    used <- str_1st(model_traits)

    use_parameters <- c(
        "beta",
        if ("l" %in% used) "LP",
        if ("d" %in% used) "DP",
        if ("t" %in% used) "RP",
        if ("l" %in% used && LP_dist == "gamma") "LP_shape",
        if ("d" %in% used && DP_dist == "gamma") "DP_shape",
        if ("t" %in% used && RP_dist == "gamma") "RP_shape",
        if (group_effect > 0) "sigma"
    )
    priors[parameter %in% use_parameters, use := TRUE]


    # Weight should not be nested unless setup is for 2+ trials
    if (!str_detect(setup, "_12")) {
        params$weight_is_nested <- FALSE
    }

    # Only use traits if in used
    GE_traits <- intersect(used, str_chars(IEs))

    pwalk(expand.grid(x = used, y = used), \(x, y) {
        priors[str_ends(parameter, str_glue("_[GE]_{x}{y}")),
               use := all(c(x, y) %in% GE_traits)]
    })

    # Now handle FEs

    fe_str <- function(base) {
        # fe_str("trial") -> "trial_[ildt]"
        # fe_str("weight") -> "weight[12]_[sildt]"
        base_str <- FEs[[base]]

        if (base_str %in% c("", "none")) return("XYZ")

        used |>
            intersect(str_chars(base_str)) |>
            str_flatten() |>
            str_glue("{base}{wt}_[{fe}]",
                     fe = _,
                     wt = if (base == "weight" && params$weight_is_nested)
                         "[12]" else "")
    }

    priors[str_detect(parameter, fe_str("trial")),  use := TRUE]
    priors[str_detect(parameter, fe_str("donor")),  use := TRUE]
    priors[str_detect(parameter, fe_str("txd")),    use := TRUE]
    priors[str_detect(parameter, fe_str("weight")), use := TRUE]


    params
}

