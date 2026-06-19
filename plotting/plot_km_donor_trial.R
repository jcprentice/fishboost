{
    library(data.table)
    library(stringr)
    library(purrr)
    library(ggplot2)
    library(HDInterval)
}

plot_km_donor_trial <- function(data_list, plot_opts = NULL) {

    if (FALSE) {
        km_data  <- readRDS("datasets/fb-test/meta/km_data_ps.rds")
        i        <- 7
        data     <- copy(km_data[[i]]$data)
        params   <- km_data[[i]]$params
        opts     <- km_data[[i]]$opts
        plot_opts <- c("keep_small_groups", "extremes", "drop_donors",
                       "mean", "fb_only", "ribbon", "t1", "t2")[c(4, 6)]
        DEBUG    <- TRUE
    } else {
        data   <- copy(data_list$data)
        params <- data_list$params
        opts   <- data_list$opts
        DEBUG  <- FALSE
    }

    if ("mean" %notin% plot_opts) {
        plotops <- setdiff(plot_opts, "ribbon")
    }

    if ("fb_only" %in% plot_opts) {
        if (DEBUG) message("- Dropping all simulated values")
        data <- data[src == "fb"]
    }

    description <- params$description |>
        str_split_1(", ") |>
        str_subset("convergence|coverage|pedigree|GRM", negate = TRUE) |>
        str_flatten_comma()

    # Choose between actual Tinfs and SIRE's inferred values for FB data
    if ("use_inferred_Tinfs" %in% plot_opts && "Tinf_sire" %in% names(data)) {
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
    if ("t1" %in% plot_opts && "t2" %notin% plot_opts) {
        data <- data[trial == 1]
    } else if ("t1" %notin% plot_opts && "t2" %in% plot_opts) {
        data <- data[trial == 2]
    }

    tmax <- params$tmax

    # If Tsign is missing but Tdeath is not, then let Tsign = Tdeath, otherwise
    # set to tmax + 1 (lines should extend past the edge of the figure)
    data[is.na(Tsign), Tsign := fifelse(!is.na(Tdeath), Tdeath, tmax[trial] + 1)]
    # RP can be extended from Tsign to tmax + 1 if Tdeath is missing
    data[, RP := Tdeath - Tsign]
    data[is.na(RP), RP := tmax[trial] - Tsign]
    data[, RP := pmax(RP, 0)]

    if ("drop_donors" %in% plot_opts) {
        if (DEBUG) message("- Dropping Donors")
        data <- data[donor == 0]
    }

    sv_curve <- function(x) c(0, sort(x, na.last = TRUE))

    # Make a copy of the data and add in rows representing t=0, sort times, make
    # sure to include NA values last.
    data_t0 <- data[!is.na(sire),
                    .(survival = seq(1, 0, length.out = .N + 1),
                      Tinf  = sv_curve(Tinf),
                      Tsign = sv_curve(Tsign),
                      RP    = sv_curve(RP),
                      src   = c(first(src), src)),
                    .(id, donor, trial)] |>
        setorder(id, donor, trial)

    # Rename src
    # data_t0[, src := str_c(src, trial)]

    # Create a column to group by
    data_t0[, gp := .GRP, .(id, donor, trial)]

    # Melt so we can use facet_wrap
    data_t1 <- melt(data_t0,
                    measure.vars = c("Tinf", "Tsign", "RP"),
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

    data_t2[variable == "Tsign",
             survival := nafill(survival, type = "locf")]
    data_t2[variable == "RP" & src == "sim",
             survival := nafill(survival, type = "locf")]
    data_t2[, tmp := seq(.N), id]
    ids_to_keep <- data_t2[src == "fb" & !is.na(survival), tmp]

    data_t2 <- data_t2[tmp %in% ids_to_keep]
    data_t2[, tmp := NULL]

    # Apply means
    if ("mean" %in% plot_opts) {
        if (DEBUG) message("- Reducing to mean of simulated curves")
        data_t2a <- data_t2[, map_if(.SD, is.numeric, mean, .else = first),
                 .(trial, donor, src, variable, time),
                 .SDcols = -"id"]

        data_t2a <- data_t2[, .(survival = mean(survival),
                                hdi1 = hdi(survival)[["lower"]],
                                hdi2 = hdi(survival)[["upper"]]),
                            .(trial, donor, src, variable, time)]

        data_t2a[, id := .GRP, src]
        data_t2a[, gp := .GRP, .(id, donor, trial)]
        setcolorder(data_t2a, names(data_t2))

        data_t3 <- data_t2a
    } else {
        data_t3 <- copy(data_t2)
    }

    data_t3 <- data_t3[time <= tmax[trial]]

    if ("show_Tinfs" %notin% plot_opts) {
        data_t3 <- data_t3[variable != "Tinf"]
    }

    data_t3[, str := str_c(src, fifelse(donor == 1, "_d", "_r"))]

    scm <- rowwiseDT(
        breaks=, labels=,              values=,
        "fb_d",  "Data Seeder",        "#33A02C",
        "sim_d", "Simulation Seeder",  "#B2DF8A",
        "fb_r",  "Data Contact",       "#1F78B4",
        "sim_r", "Simulation Contact", "#A6CEE3"
    )

    slw <- if ("mean" %in% plot_opts) 0.7 else 0.2
    slt <- if ("mean" %in% plot_opts) "dashed" else "solid"

    # Plot
    ribbon <- if (all(c("mean", "ribbon") %in% plot_opts)) {
        geom_ribbon(aes(x = time, ymin = hdi1, ymax = hdi2,
                group = gp, fill = str, linewidth = src),
            data_t3[src == "sim"],
            alpha = 0.5, show.legend = FALSE)
    }

    plt <- ggplot() +
        ribbon +
        geom_line(aes(x = time,
                      y = survival,
                      group = gp,
                      colour = str,
                      linewidth = src,
                      linetype = src),
                  data_t3) +
        scale_fill_manual(breaks = scm$breaks,
                          labels = scm$labels,
                          values = scm$values) +
        scale_colour_manual(breaks = scm$breaks,
                            labels = scm$labels,
                            values = scm$values) +
        scale_linewidth_manual(breaks = c("fb", "sim"),
                               values = c(0.5, slw),
                               guide = "none") +
        scale_linetype_manual(breaks = c("fb", "sim"),
                              values = c("solid", slt),
                              guide = "none") +
        lims(y = 0:1) +
        # coord_cartesian(expand = FALSE) +
        labs(colour = "Source",
             #linewidth = "Source",
             x = "Time (days)",
             y = "Proportion",
             title = str_glue("KM by seeder status, {opts$post}"),
             subtitle = str_glue("{params$dataset}/{params$label}: {description}")) +
        facet_grid(cols = vars(variable),
                   rows = vars(trial),
                   scales = "free_x",
                   labeller = labeller(
                       variable = c(Tinf  = "Time to infection",
                                    Tsign = "Time to visual signs",
                                    RP    = "Time from visual signs to death"),
                       trial = c("1" = "Trial 1",
                                 "2" = "Trial 2"))) +
        theme_bw() +
        theme(panel.background = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              legend.position = "bottom")
    plt

    list(plt = plt, data = data_t3)

}

