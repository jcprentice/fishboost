{
    library(data.table)
    library(stringr)
    library(purrr)
    library(ggplot2)
    library(HDInterval)
}

plot_km_sire_trial <- function(data_list, plotopts = NULL) {

    if (FALSE) {
        km_data  <- readRDS("datasets/fb-test/meta/km_data_ps.rds")
        i        <- 7
        data     <- copy(km_data[[i]]$data)
        params   <- km_data[[i]]$params
        opts     <- km_data[[i]]$opts
        plotopts <- c("keep_small_groups", "extremes", "drop_donors",
                      "mean", "fb_only", "ribbon", "t1", "t2")[c(4, 6)]
        DEBUG    <- TRUE
    } else {
        data   <- copy(data_list$data)
        params <- data_list$params
        opts   <- data_list$opts
        DEBUG  <- FALSE
    }

    if ("mean" %notin% plotopts) {
        plotops <- setdiff(plotopts, "ribbon")
    }

    if ("fb_only" %in% plotopts) {
        if (DEBUG) message("- Dropping all simulated values")
        data <- data[src == "fb"]
    }

    # Donor families are any families
    data[, donor := fifelse(any(donor == 1L), 1L, 0L), .(sire, trial)]

    description <- params$description |>
        str_split_1(", ") |>
        str_subset("convergence|coverage|pedigree|GRM", negate = TRUE) |>
        str_replace_all(c("inf_model 1" = "inf: I = D",
                          "inf_model 2" = "inf: I = 0.1*D",
                          "inf_model 3" = "inf: Don = 0.1*Rec",
                          "inf_model 4" = "inf: Don = r*Rec")) |>
        str_flatten_comma()

    # Choose between actual Tinfs and SIRE's inferred values for FB data
    if ("use_inferred_Tinfs" %in% plotopts && "Tinf_sire" %in% names(data)) {
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
    if ("t1" %in% plotopts && "t2" %notin% plotopts) {
        data <- data[trial == 1]
    } else if ("t1" %notin% plotopts && "t2" %in% plotopts) {
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

    if ("drop_donors" %in% plotopts) {
        if (DEBUG) message("- Dropping Donors")
        data <- data[donor == 0]
    }

    sv_curve <- function(x) c(0, sort(x, na.last = TRUE))

    # Make a copy of the data and add in rows representing t=0, sort times, make
    # sure to include NA values last
    data_t0 <- data[!is.na(sire),
                    .(donor = fifelse(any(donor == 1), 1, 0),
                      survival = seq(1, 0, length.out = .N + 1),
                      Tinf  = sv_curve(Tinf),
                      Tsign = sv_curve(Tsign),
                      RP    = sv_curve(RP),
                      src   = c(first(src), src)),
                    .(id, sire, trial)] |>
        setorder(id, sire, trial)

    # Create a column to group by
    data_t0[, gp := .GRP, .(id, sire, trial, donor)]

    if ("keep_small_groups" %notin% plotopts) {
        if (DEBUG) message("- Dropping small groups")
        small_groups <- data_t0[, .(N = .N), gp][N < 10, gp]
        data_t0 <- data_t0[gp %notin% small_groups]
        data_t0[, gp := .GRP, .(id, sire, trial, donor)]
    }

    if ("extremes" %in% plotopts) {
        if (DEBUG) message("- Keeping only extremes")
        foo <- data_t0[id == last(id), .(sire, trial, donor, Tsign)]
        # Filter out any sires with fewer than 5 non-NA values
        foo[, p := .N - sum(is.na(Tsign)), .(sire, trial)]
        foo <- foo[p > 5]
        foo1 <- foo[, .(mu = mean(Tsign, na.rm = TRUE)),
                    .(sire, trial, donor)] |>
            setorder(trial, sire, mu)
        # foo1[, donor := fifelse(sire %in% donor_sires, 1L, 0L)]
        ids <- foo1[, .(sire = sire[c(which.min(mu), which.max(mu))],
                        es = c("lo", "hi")),
                    .(trial, donor)]
        data_t0 <- merge(data_t0, ids, by = c("sire", "trial", "donor")) |>
            setcolorder("donor", after = "trial")
    } else {
        data_t0[, es := "X"]
    }

    # Melt so we can use facet_wrap
    data_t1 <- melt(data_t0,
                    measure.vars = c("Tinf", "Tsign", "RP"),
                    value.name = "time") |>
        setcolorder("survival", after = "time")


    # Use approxfun to convert to t = 0,1,...,tmax
    data_t2 <- data_t1[, {
        times <- seq(0, tmax[trial[[1]]])
        list(es = first(es),
             time = times,
             survival = approx(time, survival, times,
                               method = "constant",
                               ties = list("ordered", max))$y)
    },
    .(id, sire, trial, donor, src, gp, variable)]

    data_t2[variable == "Tsign" | (variable == "RP" & src == "sim"),
            survival := nafill(survival, type = "locf")]
    data_t2[, tmp := seq(.N), id]
    ids_to_keep <- data_t2[src == "fb" & !is.na(survival), tmp]

    data_t2 <- data_t2[tmp %in% ids_to_keep, .SD]
    data_t2[, tmp := NULL]


    # Apply means and HDI
    if ("mean" %in% plotopts) {
        if (DEBUG) message("- Reducing to mean of simulated curves")
        data_t2a <- if ("extremes" %in% plotopts) {
            data_t2[, .(es = first(es),
                        survival = mean(survival),
                        hdi1 = hdi(survival)[["lower"]],
                        hdi2 = hdi(survival)[["upper"]]),
                    .(sire, trial, donor, src, variable, time)]
        } else {
            rbind(
                data_t2[src == "sim",
                        .(sire = 0,
                          es = first(es),
                          survival = mean(survival),
                          hdi1 = hdi(survival)[["lower"]],
                          hdi2 = hdi(survival)[["upper"]]),
                        .(trial, donor, src, variable, time)],
                data_t2[src == "fb",
                        .(gp = first(gp),
                          es = first(es),
                          survival = mean(survival)),
                        .(sire, trial, donor, src, variable, time)],
                fill = TRUE) |>
                setorder(-src, trial, sire, donor, variable, time)
        }
        data_t2a[, id := .GRP, src]
        data_t2a[, gp := .GRP, .(id, sire, trial, donor)]
        setcolorder(data_t2a, names(data_t2))

        data_t3 <- data_t2a
    } else {
        data_t3 <- copy(data_t2)
    }

    data_t3 <- data_t3[time <= tmax[trial]]

    if ("show_Tinfs" %notin% plotopts) {
        data_t3 <- data_t3[variable != "Tinf"]
    }

    data_t3[, str := str_c(src, es, fifelse(donor == 1, "d", "r"), sep = "_") |>
                str_remove("X_")]

    scm <- rowwiseDT(
        breaks=,    labels=,                   values=,
        "fb_hi_d",  "Data Seeder high",        "#33A02C",
        "sim_hi_d", "Simulation Seeder high",  "#B2DF8A",
        "fb_lo_d",  "Data Seeder low",         "#E31A1C",
        "sim_lo_d", "Simulation Seeder low",   "#FB9A99",
        "fb_hi_r",  "Data Contact high",       "#1F78B4",
        "sim_hi_r", "Simulation Contact high", "#A6CEE3",
        "fb_lo_r",  "Data Contact low",        "#FF7F00",
        "sim_lo_r", "Simulation Contact low",  "#FDBF6F",
        "fb_d",     "Data Seeder",             "#33A02C",
        "sim_d",    "Simulation Seeder",       "#B2DF8A",
        "fb_r",     "Data Contact",            "#1F78B4",
        "sim_r",    "Simulation Contact",      "#A6CEE3"
    )

    if ("extremes" %notin% plotopts) {
        scm[, labels := str_remove(labels, " high| low")]
    }

    slw <- if ("mean" %in% plotopts) 0.7 else 0.2
    slt <- if ("mean" %in% plotopts) "dashed" else "solid"

    # Plot
    plt <- ggplot() +
        {if (all(c("mean", "ribbon") %in% plotopts)) {
            geom_ribbon(aes(x = time, ymin = hdi1, ymax = hdi2,
                            group = gp, fill = str, linewidth = src),
                        data_t3[src == "sim"],
                        alpha = 0.5, show.legend = FALSE)}} +
        geom_line(aes(x = time, y = survival,
                      group = gp, colour = str,
                      linewidth = src, linetype = src),
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
        lims(y = c(0, 1)) +
        # coord_cartesian(expand = FALSE) +
        labs(colour = "Source",
             #linewidth = "Source",
             x = "Time (days)",
             y = "Proportion",
             title = str_glue("KM by sire family, {opts$post}"),
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

    list(plt = plt, data = data_t2)
}
