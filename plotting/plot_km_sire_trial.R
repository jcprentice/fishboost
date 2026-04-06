{
    library(data.table)
    library(stringr)
    library(purrr)
    library(ggplot2)
}

plot_km_sire_trial <- function(data_list, plotopts = NULL) {

    if (FALSE) {
        i <- 1
        data <- copy(km_data[[i]]$data)
        params <- km_data[[i]]$params
        opts <- km_data[[i]]$opts
        plotopts <- c("drop_small_groups", "extreme_sires", "drop_donors",
                      "mean", "t1", "t2")[c(1, 4)]
    } else {
        data   <- copy(data_list$data)
        params <- data_list$params
        opts   <- data_list$opts
    }

    # Get trials (need as character for indexing donor_sires)
    # trials <- data[, unique(trial)] |> discard(is.na) |> as.character()

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
    # set to tmax + 1
    data[is.na(Tsym), Tsym := fifelse(!is.na(Tdeath), Tdeath, tmax[trial] + 1)]
    # RP can be extended from Tsym to tmax + 1 if Tdeath is missing
    data[, RP := Tdeath - Tsym]
    data[is.na(RP), RP := tmax[trial] - Tsym]

    if ("drop_donors" %in% plotopts) {
        data <- data[donor == 0]
    }

    sv_curve <- function(x) c(0, sort(x, na.last = TRUE))

    # Make a copy of the data and add in rows representing t=0, sort times, make
    # sure to include NA values last
    data_t0 <- data[!is.na(sire),
                    .(donor = fifelse(any(donor == 1), 1, 0),
                      Tinf = sv_curve(Tinf),
                      Tsym = sv_curve(Tsym),
                      RP   = sv_curve(RP),
                      src = src[c(1, 1:.N)]),
                    .(id, sire, trial)] |>
        setorder(id, sire, trial)

    # Rename src
    # data_t0[, src := str_c(src, trial)]

    # Create a column to group by
    data_t0[, gp := .GRP, .(id, sire, trial)]

    # Survival curves by id and sire
    data_t0[, survival := seq(1, 0, length.out = .N), gp]

    if ("drop_small_groups" %in% plotopts) {
        data_t0 <- data_t0[, gpsize := .N, gp][gpsize >= 10]
        data_t0[, gpsize := NULL]
    }

    if ("extreme_sires" %in% plotopts) {
       foo <- data_t0[id == last(id), .(sire, trial, Tsym)]
       # Filter out any sires with fewer than 5 non-NA values
       foo[, p := .N - sum(is.na(Tsym)), .(sire, trial)]
       foo <- foo[p > 5]
       foo1 <- foo[, .(mu = mean(Tsym, na.rm = TRUE)), .(sire, trial)] |>
           setorder(trial, sire, mu)
       foo1[, donor := fifelse(sire %in% unlist(donor_sires), 1L, 0L)]
       ids <- foo1[, .(sire = sire[c(which.min(mu), which.max(mu))],
                       str = str_c(c("lo", "hi"), "_", fifelse(donor == 1, "d", "r"))),
                   .(trial, donor)]
       data_t0 <- merge(data_t0, ids, by = c("sire", "trial", "donor")) |>
           setcolorder("donor", after = "trial")
    }

    # Melt so we can use facet_wrap
    data_t1 <- melt(data_t0,
                    measure.vars = c("Tinf", "Tsym", "RP"),
                    value.name = "time")
    setcolorder(data_t1, "survival", after = "time")


    # Use approxfun to convert to t = 0,1,...,tmax
    data_t2 <- data_t1[, {
        times <- seq(0, tmax[trial[[1]]])
        list(time = times,
             survival = approx(time, survival, times,
                               method = "constant",
                               ties = list("ordered", max))$y)
    },
    .(id, sire, trial, donor, src, gp, variable)]

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
        data_t2a <- data_t2[, map(.SD, mean),
                 .(sire, trial, donor, src, variable, time),
                 .SDcols = -"id"]
        data_t2a[, id := .GRP, src]
        data_t2a[, gp := .GRP, .(id, sire, trial)]
        setcolorder(data_t2a, names(data_t2))

        data_t2 <- data_t2a
    }

    data_t2 <- data_t2[time <= tmax[trial]]

    # Drop the really high values
    # data_t2[, time := fifelse(time > c(104, 160)[trial], NA, time)]
    # data_t2[(trial == 1 & time >= 104) | (trial == 2 & time >= 160), time := NA]
    # maxvals <- data[, ceiling(max(c(Tsym, Tdeath), na.rm = TRUE)), trial][, V1]
    # maxvals <- c(104, 160)
    # data_t2[, time := fifelse(time >= maxvals[trial], maxvals[trial], time)]

    # we might want to get a subset of the best and worst sires in each trial
    # if ("extremes" %in% plotopts) {
    #     sires1 <- (data_t2[str_starts(src, "fb1") & variable == "Tsym" & !is.na(time),
    #                        .(survival = survival[.N]), sire]
    #                [, .(sire_min = sire[which.min(survival)],
    #                     sire_max = sire[which.max(survival)])]
    #                [, c(sire_min, sire_max) |> sort()])
    #
    #     sires2 <- (data_t2[str_starts(src, "fb2") & variable == "Tsym" & !is.na(time),
    #                        .(survival = survival[.N]), sire]
    #                [, .(sire_min = sire[which.min(survival)],
    #                     sire_max = sire[which.max(survival)])]
    #                [, c(sire_min, sire_max) |> sort()])
    #
    #     data_t2 <- data_t2[(trial == 1 & sire %in% sires1) | (trial == 2 & sire %in% sires2)]
    # }

    if ("show_Tinfs" %notin% plotopts) {
        data_t2 <- data_t2[variable != "Tinf"]
    }

    if ("extreme_sires" %in% plotopts) {
        data_t2[, str := str_c(src, "_", str)]
    } else {
        data_t2[, str := src]
    }

    breaks <- c("fb",  "fb_lo_d",  "fb_hi_d",  "fb_lo_r",  "fb_hi_r",
                "sim", "sim_lo_d", "sim_hi_d", "sim_lo_r", "sim_hi_r")

    labels <- c("Data",
                "Data Seeder low",  "Data Seeder high",
                "Data Contact low", "Data Contact high",
                "Simulation",
                "Simulation Seeder low",  "Simulation Seeder high",
                "Simulation Contact low", "Simulation Contact high")

    #             Regular    Low donor  High donor Low recip  High recip
    col_vals <- c("#1F78B4", "#FF7F00", "#1F78B4", "#E31A1C", "#33A02C", # Data
                  "#A6CEE3", "#FDBF6F", "#A6CEE3", "#FB9A99", "#B2DF8A") # Simulation

    if ("mean" %in% plotopts) {
        col_vals <- c("#1F78B4", "#FF7F00", "#1F78B4", "#E31A1C", "#33A02C", # Data
                      "#E31A1C", "#FDBF6F", "#A6CEE3", "#FB9A99", "#B2DF8A") # Simulation
    }


    # Plot
    plt <- ggplot(data_t2) +
        # geom_step
        geom_line(aes(x = time,
                      y = survival,
                      group = gp,
                      colour = str,
                      linewidth = src)) +
        scale_colour_manual(breaks = breaks,
                            labels = labels,
                            values = col_vals) +
        scale_linewidth_manual(breaks = c("fb", "sim"),
                               values = c(0.5, 0.2),
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

    # gfx_dir <- str_glue("{dataset}/gfx")

    # ggsave(str_glue("{gfx_dir}/{dataset}-s{scen}-KM_grid-trials.pdf"),
    #        plot = plt_trials, width = 12, height = 8)
    # ggsave(str_glue("{gfx_dir}/{dataset}-s{scen}-KM_grid-trials.png"),
    #        plot = plt_trials, width = 12, height = 8, dpi = 600)
    #
    # plt2 <- plt +
    #     facet_wrap(. ~ variable,
    #                scales = "free_x",
    #                labeller = labeller(
    #                    variable = c(Tinf = "Proportion of family uninfected vs time",
    #                                 Tsym = "Proportion of family with no symptoms vs time",
    #                                 RP   = "Proportion of family surviving vs time")))
    #
    # plt2

    # ggsave(str_glue("{gfx_dir}/{dataset}-s{scen}-KM_grid.pdf"),
    #        plot = plt2, width = 12, height = 6)
    # ggsave(str_glue("{gfx_dir}/{dataset}-s{scen}-KM_grid.png"),
    #        plot = plt2, width = 12, height = 6, dpi = 600)

    list(plt = plt, data = data_t2)
}
