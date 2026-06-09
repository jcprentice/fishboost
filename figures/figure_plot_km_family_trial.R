{
    library(data.table)
    library(stringr)
    library(purrr)
    library(ggplot2)
    library(HDInterval)
}

figure_plot_km_family_trial <- function() {

    {
        km_data  <- readRDS("datasets/fb-test/meta/km_data_ps.rds")
        i        <- 2
        data     <- copy(km_data[[i]]$data)
        params   <- km_data[[i]]$params
        opts     <- km_data[[i]]$opts
        plotopts <- c("keep_small_groups", "extremes", "drop_donors", "mean", "ribbon")
        plotopts <- c("ribbon")
    }

    # Set families
    data[, family := .GRP, .(sire, dam, trial)] |>
        setcolorder("family", after = "id")

    # Choose between actual Tinfs and inferred values for FB data
    if ("use_inferred_Tinfs" %in% plotopts && "Tinf_sire" %in% names(data)) {
        message("- Using sire Tinfs")
        # data[src == "fb", Tinf := Tinf_sire]
        median_Tinfs <- get_median_Tinfs(params)

        if (length(median_Tinfs) == data[src == "fb", .N]) {
            data[src == "fb", Tinf := median_Tinfs]
        }
    }

    # Drop small families unless specified otherwise
    if ("keep_small_groups" %notin% plotopts) {
        message("- Dropping small groups")
        small_groups <- data[id == last(id), .N, family][N < 10, family]
        data <- data[family %notin% small_groups]
    } else {
        message("- Keeping small groups")
    }

    # Donor families are any families
    data[, donor := fifelse(any(donor == 1L), 1L, 0L), .(family, trial)]

    # Filter by trial
    if ("t1" %in% plotopts & "t2" %notin% plotopts) {
        message(" - Keeping Trial 1 only")
        data <- data[trial == 1]
    } else if ("t1" %notin% plotopts & "t2" %in% plotopts) {
        message(" - Keeping Trial 2 only")
        data <- data[trial == 2]
    }

    tmax <- params$tmax

    # If Tsign is missing but Tdeath is not, then let Tsign = Tdeath, otherwise
    # set to tmax + 10 (lines should extend past the edge of the figure)
    data[is.na(Tsign), Tsign := fifelse(!is.na(Tdeath), Tdeath, tmax[trial] + 1)]
    # RP can be extended from Tsign to tmax + 1 if Tdeath is missing
    data[, RP := Tdeath - Tsign]
    data[is.na(RP), RP := tmax[trial] - Tsign]
    data[, RP := pmax(RP, 0)]

    if ("drop_donors" %in% plotopts) {
        message("- Dropping Donors")
        data <- data[donor == 0]
    }

    sv_curve <- function(x) c(0, sort(x, na.last = TRUE))

    # Make a copy of the data and add in rows representing t=0, sort times, make
    # sure to include NA values last
    data_t0 <- data[!is.na(family),
                    .(donor = fifelse(any(donor == 1), 1, 0),
                      survival = seq(1, 0, length.out = .N + 1),
                      Tinf  = sv_curve(Tinf),
                      Tsign = sv_curve(Tsign),
                      RP    = sv_curve(RP)),
                    .(id, src, family, trial)] |>
        setorder(id, family, trial)

    # Create a column to group by
    data_t0[, gp := .GRP, .(id, family, trial, donor)]

    # Set extreme families
    foo <- data_t0[id == last(id), .(family, trial, donor, Tsign)]
    # Filter out any sires with fewer than 5 non-NA values
    foo[, p := .N - sum(is.na(Tsign)), .(family, trial)]
    foo <- foo[p > 5]
    foo1 <- foo[, .(mu = mean(Tsign, na.rm = TRUE)),
                .(family, trial, donor)] |>
        setorder(trial, family, mu)
    # foo1[, donor := fifelse(family %in% donor_sires, 1L, 0L)]
    ids <- foo1[, .(family = family[c(which.min(mu), which.max(mu))],
                    es = c("lo", "hi")),
                .(trial, donor)]
    data_t0[, es := fcase(
        id %in% ids[es == "lo", family], "lo",
        id %in% ids[es == "hi", family], "hi",
        default = "X"
    )]

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
    .(id, family, trial, donor, src, gp, variable)]

    data_t2[variable == "Tsign" | (variable == "RP" & src == "sim"),
            survival := nafill(survival, type = "locf")]
    data_t2[, tmp := seq(.N), id]
    ids_to_keep <- data_t2[src == "fb" & !is.na(survival), tmp]

    data_t2 <- data_t2[tmp %in% ids_to_keep, .SD]
    data_t2[, tmp := NULL]


    # Apply means and HDI
    if ("mean" %in% plotopts) {
        message("- Reducing to mean of simulated curves")
        data_t2a <- if ("extremes" %in% plotopts) {
            data_t2[, .(es = first(es),
                        survival = mean(survival),
                        hdi1 = hdi(survival)[["lower"]],
                        hdi2 = hdi(survival)[["upper"]]),
                    .(family, trial, donor, src, variable, time)]
        } else {
            rbind(
                data_t2[src == "sim",
                        .(family = 0,
                          es = first(es),
                          survival = mean(survival),
                          hdi1 = hdi(survival)[["lower"]],
                          hdi2 = hdi(survival)[["upper"]]),
                        .(trial, donor, src, variable, time)],
                data_t2[src == "fb",
                        .(gp = first(gp),
                          es = first(es),
                          survival = mean(survival)),
                        .(family, trial, donor, src, variable, time)],
                fill = TRUE) |>
                setorder(-src, trial, family, donor, variable, time)
        }
        data_t2a[, id := .GRP, src]
        data_t2a[, gp := .GRP, .(id, family, trial, donor)]
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

    # scm <- rowwiseDT(
    #     breaks=, labels=,          values=,
    #     "fb_d",  "Seeder family",  "#33A02C",
    #     "fb_r",  "Contact family", "#1F78B4"
    # )

    if ("extremes" %notin% plotopts) {
        scm[, labels := str_remove(labels, " high| low")]
    }

    slw <- if ("mean" %in% plotopts) 0.7 else 0.2
    slt <- if ("mean" %in% plotopts) "dashed" else "solid"

    description <- params$description |>
        str_split_1(", ") |>
        str_subset("convergence|coverage|pedigree|GRM", negate = TRUE) |>
        str_flatten_comma()


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
             title = str_glue("KM by full-sib family, {opts$post}"),
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
