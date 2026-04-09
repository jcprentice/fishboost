{
    library(data.table)
    library(stringr)
    library(purrr)
    library(ggplot2)
}

plot_km_donor_trial <- function(data_list, plotopts = NULL) {

    if (FALSE) {
        km_data  <- readRDS("datasets/fb-test/meta/km_data_ps.rds")
        i        <- 1
        data     <- copy(km_data[[i]]$data)
        params   <- km_data[[i]]$params
        opts     <- km_data[[i]]$opts
        plotopts <- c("drop_small_groups", "extreme_sires", "drop_donors",
                      "mean", "fb_only", "t1", "t2")[c(2, 4)]
        DEBUG    <- TRUE
    } else {
        data   <- copy(data_list$data)
        params <- data_list$params
        opts   <- data_list$opts
        DEBUG  <- FALSE
    }

    if ("fb_only" %in% plotopts) {
        if (DEBUG) message("- Dropping all simulated values")
        data <- data[src == "fb"]
    }

    # Families who have any donors
    donor_sires <- data[donor == 1, sort(unique(sire)), trial] |>
        split(by = "trial") |> map("V1")

    description <- params$description |>
        str_split_1(", ") |>
        str_subset("convergence|coverage|pedigree|GRM", negate = TRUE) |>
        str_replace_all(c("inf_model 1" = "inf: I = D",
                          "inf_model 2" = "inf: I = 0.1*D",
                          "inf_model 3" = "inf: Don = 0.1*Rec",
                          "inf_model 4" = "inf: Don = r*Rec")) |>
        str_flatten_comma()

    # Choose between actual Tinfs and SIRE's inferred values for FB data
    if ("use_sire_Tinfs" %in% plotopts && "Tinf_sire" %in% names(data)) {
        if (DEBUG) message("- Using sire Tinfs")
        # data[src == "fb", Tinf := Tinf_sire]
        Tinf_names <- str_c("Tinf_sire_", seq_len(opts$n_plots))
        fb_Tinfs <- with(params, get_Tinfs(dataset, scenario, replicate, opts$n_plots)) |>
            as.data.table() |>
            setnames(Tinf_names)

        fb_Tinfs <- cbind(data[src == "fb", -c("Tinf")],
                          fb_Tinfs[55:1829])
    }

    # Filter by trial
    if ("t1" %in% plotopts) data <- data[trial == 1]
    if ("t2" %in% plotopts) data <- data[trial == 2]

    tmax <- params$tmax

    # If Tsym is missing but Tdeath is not, then let Tsym = Tdeath, otherwise
    # set to tmax + 1 (lines should extend past the edge of the figure)
    data[is.na(Tsym), Tsym := fifelse(!is.na(Tdeath), Tdeath, tmax[trial] + 1)]
    # RP can be extended from Tsym to tmax + 1 if Tdeath is missing
    data[, RP := Tdeath - Tsym]
    data[is.na(RP), RP := tmax[trial] - Tsym]
    data[, RP := pmax(RP, 0)]

    if ("drop_donors" %in% plotopts) {
        if (DEBUG) message("- Dropping Donors")
        data <- data[donor == 0]
    }

    sv_curve <- function(x) c(0, sort(x, na.last = TRUE))

    # Make a copy of the data and add in rows representing t=0, sort times, make
    # sure to include NA values last.
    data_t0 <- data[, .(survival = seq(1, 0, length.out = .N + 1),
                        Tinf = sv_curve(Tinf),
                        Tsym = sv_curve(Tsym),
                        RP   = sv_curve(RP),
                        src = c(first(src), src)),
                    .(id, donor, trial)] |>
        setorder(id, donor, trial)

    # Rename src
    # data_t0[, src := str_c(src, trial)]

    # Create a column to group by
    data_t0[, gp := .GRP, .(id, donor, trial)]

    if ("drop_small_groups" %in% plotopts) {
        if (DEBUG) message("- Dropping small groups")
        small_groups <- data_t0[, .(N = .N), gp][N < 10, gp]
        data_t0 <- data_t0[gp %notin% small_groups]
        data_t0[, gp := .GRP, .(id, sire, trial)]
    }

    # Melt so we can use facet_wrap
    data_t1 <- melt(data_t0,
                    measure.vars = c("Tinf", "Tsym", "RP"),
                    value.name = "time") |>
        setcolorder("survival", after = "time")


    # Use approxfun to convert to t = 0,1,...,tmax
    data_t2 <- data_t1[, {
        times <- seq(0, tmax[trial[[1]]])
        list(time = times,
             survival = approx(time, survival, times,
                               method = "constant",
                               ties = list("ordered", max))$y)
    },
    .(id, trial, donor, src, gp, variable)]

    data_t2[variable == "Tsym",
             survival := nafill(survival, type = "locf")]
    data_t2[variable == "RP" & src == "sim",
             survival := nafill(survival, type = "locf")]
    data_t2[, tmp := seq(.N), id]
    ids_to_keep <- data_t2[src == "fb" & !is.na(survival), tmp]

    data_t2 <- data_t2[tmp %in% ids_to_keep, .SD]
    data_t2[, tmp := NULL]

    # Apply means
    if ("mean" %in% plotopts) {
        if (DEBUG) message("- Reducing to mean of simulated curves")
        data_t2a <- data_t2[, map_if(.SD, is.numeric, mean, .else = first),
                 .(trial, donor, src, variable, time),
                 .SDcols = -"id"]
        data_t2a[, id := .GRP, src]
        data_t2a[, gp := .GRP, .(id, donor, trial)]
        setcolorder(data_t2a, names(data_t2))

        data_t2 <- data_t2a
    }

    data_t2 <- data_t2[time <= tmax[trial]]

    if ("show_Tinfs" %notin% plotopts) {
        data_t2 <- data_t2[variable != "Tinf"]
    }

    data_t2[, str := str_c(src, fifelse(donor == 1, "_d", "_r"))]

    scm <- rowwiseDT(
        breaks=,    labels=,                   values=,
        "fb_hi_d",  "Data Seeder high",        "#33A02C",
        "sim_hi_d", "Simulation Seeder high",  "#B2DF8A",
        "fb_lo_d",  "Data Seeder low",         "#E31A1C",
        "fb_hi_r",  "Data Contact high",       "#1F78B4",
        "sim_lo_d", "Simulation Seeder low",   "#FB9A99",
        "sim_hi_r", "Simulation Contact high", "#A6CEE3",
        "fb_lo_r",  "Data Contact low",        "#FF7F00",
        "sim_lo_r", "Simulation Contact low",  "#FDBF6F"
    )

    if ("extreme_sires" %notin% plotopts) {
        scm[, labels := str_remove(labels, " high| low")]
    }

    slw <- if ("mean" %in% plotopts) 0.7 else 0.2

    # Plot
    plt <- ggplot(data_t2) +
        # geom_step
        geom_line(aes(x = time,
                      y = survival,
                      group = gp,
                      colour = str,
                      linewidth = src)) +
        scale_colour_manual(breaks = scm$breaks,
                            labels = scm$labels,
                            values = scm$values) +
        scale_linewidth_manual(breaks = c("fb", "sim"),
                               values = c(0.5, slw),
                               guide = "none") +
        lims(y = 0:1) +
        # coord_cartesian(expand = FALSE) +
        labs(colour = "Source",
             #linewidth = "Source",
             x = "Time (days)",
             y = "Proportion",
             title = str_glue("KM by family (sires), {opts$post}"),
             subtitle = str_glue("{params$dataset}/{params$label}: {description}")) +
        facet_grid(cols = vars(variable),
                   rows = vars(trial),
                   scales = "free_x",
                   labeller = labeller(
                       variable = c(Tinf = "Proportion of family uninfected vs time",
                                    Tsym = "Proportion of family with no symptoms vs time",
                                    RP   = "Proportion of family surviving vs time"),
                       trial = c("1" = "Trial 1",
                                 "2" = "Trial 2"))) +
        theme_bw() +
        theme(panel.background = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              legend.position = "bottom")
    plt

    list(plt = plt, data = data_t2)

}

