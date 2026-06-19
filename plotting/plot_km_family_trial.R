{
    library(data.table)
    library(stringr)
    library(purrr)
    library(ggplot2)
    library(HDInterval)

    source("get_Tinfs.R")
}

plot_km_family_trial <- function(data_list, plot_opts = NULL) {

    if (FALSE) {
        km_data  <- readRDS("datasets/fb-test/meta/km_data_ps.rds")
        i        <- 2
        data     <- copy(km_data[[i]]$data)
        params   <- km_data[[i]]$params
        opts     <- km_data[[i]]$opts
        plot_opts <- c("keep_small_groups", "extremes", "drop_donors",
                       "mean", "fb_only", "ribbon", "t1", "t2",
                       "show_Tinfs", "use_inferred_Tinfs")[5]
        DEBUG    <- TRUE
    } else {
        data   <- copy(data_list$data)
        params <- data_list$params
        opts   <- data_list$opts
        DEBUG  <- FALSE
    }

    if ("fb_only" %in% plot_opts) {
        if (DEBUG) message("- Dropping all simulated values")
        data <- data[src == "fb"]
    }

    # Set families
    data[, family := .GRP, .(sire, dam, trial)] |>
        setcolorder("family", after = "id")

    # Choose between actual Tinfs and inferred values for FB data
    if ("use_inferred_Tinfs" %in% plot_opts) {
        if (DEBUG) message("- Using inferred Tinfs")
        # data[src == "fb", Tinf := Tinf_sire]
        median_Tinfs <- get_median_Tinfs(params)

        if (length(median_Tinfs) == data[src == "fb", .N]) {
            data[src == "fb", Tinf := median_Tinfs]
        }
    }

    # Drop small families unless specified otherwise
    if ("keep_small_groups" %notin% plot_opts) {
        small_groups <- data[id == last(id), .N, family][N < 10, family]
        data <- data[family %notin% small_groups]
    } else {
        if (DEBUG) message("- Keeping small groups")
    }

    # Filter by trial
    if ("t1" %in% plot_opts && "t2" %notin% plot_opts) {
        if (DEBUG) message(" - Keeping Trial 1 only")
        data <- data[trial == 1]
    } else if ("t1" %notin% plot_opts && "t2" %in% plot_opts) {
        if (DEBUG) message(" - Keeping Trial 2 only")
        data <- data[trial == 2]
    }

    # If Tsign is missing but Tdeath is not, then let Tsign = Tdeath, otherwise
    # set to tmax + 10 (lines should extend past the edge of the figure)
    # RP can be extended from Tsign to tmax + 1 if Tdeath is missing
    {
        tmax <- params$tmax
        data[is.na(Tsign), Tsign := fifelse(!is.na(Tdeath), Tdeath, tmax[trial] + 1)]
        data[, RP := Tdeath - Tsign]
        data[is.na(RP), RP := tmax[trial] - Tsign]
        data[, RP := pmax(RP, 0)]
    }

    # Donor families are any families
    data[, donor := fifelse(any(donor == 1L), 1L, 0L), .(family, trial)]

    if ("drop_donors" %in% plot_opts) {
        if (DEBUG) message("- Dropping Donors")
        data <- data[donor == 0]
    }

    # Make a copy of the data and add in rows representing t=0, sort times, make
    # sure to include NA values last
    sv_curve <- function(x) c(0, sort(x, na.last = TRUE))

    data_sv <- data[!is.na(family),
                    .(survival = seq(100, 0, length.out = .N + 1),
                      Tinf  = sv_curve(Tinf),
                      Tsign = sv_curve(Tsign),
                      RP    = sv_curve(RP)),
                    .(id, src, family, trial, donor)] |>
        setorder(id, family, trial)

    # Create a column to group by
    data_sv[, gp := .GRP, .(id, family, trial, donor)]

    if ("extremes" %in% plot_opts) {
        if (DEBUG) message("- Keeping only extremes")
        foo <- data_sv[id == last(id), .(family, trial, donor, Tsign)]
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
        data_sv <- merge(data_sv, ids, by = c("family", "trial", "donor")) |>
            setcolorder("donor", after = "trial")
    } else {
        data_sv[, es := "X"]
    }

    if ("show_Tinfs" %notin% plot_opts) {
        data_sv[, Tinf := NULL]
    }

    # Melt so we can use facet_wrap
    data_svm <- melt(data_sv,
                     measure.vars = patterns("Tinf|Tsign|RP"),
                     value.name = "time") |>
        setcolorder("survival", after = "time") |>
        # Use approxfun to convert to t = 0,1,...,tmax
        _[, {
            times <- seq(0, tmax[trial[[1]]])
        list(es = first(es),
             time = times,
             survival = approx(time, survival, times,
                               method = "constant",
                               ties = list("ordered", max))$y)
    },
    .(id, family, trial, donor, src, gp, variable)]

    data_svm[, survival := nafill(survival, type = "locf")]

    data_svm[, str := str_c(src, es, fifelse(donor == 1, "d", "r"), sep = "_") |>
                 str_remove("X_")]


    # Get simulation means and HDI
    data_sim <- data_svm[src == "sim",
                         .(gp = first(gp),
                           survival = mean(survival),
                           hdi1 = hdi(survival)[["lower"]],
                           hdi2 = hdi(survival)[["upper"]],
                           str = first(str)),
                         .(trial, donor, src, es, variable, time)] |>
        setcolorder(names(data_svm), skip_absent = TRUE)


    # Plot

    scm <- rowwiseDT(
        breaks=,    labels=,                   values=,
        "fb_d",     "Data Seeder",             "#33A02C",
        "sim_d",    "Simulation Seeder",       "#B2DF8A",
        "fb_r",     "Data Contact",            "#1F78B4",
        "sim_r",    "Simulation Contact",      "#A6CEE3",
        "fb_hi_d",  "Data Seeder high",        "#33A02C",
        "sim_hi_d", "Simulation Seeder high",  "#B2DF8A",
        "fb_lo_d",  "Data Seeder low",         "#E31A1C",
        "sim_lo_d", "Simulation Seeder low",   "#FB9A99",
        "fb_hi_r",  "Data Contact high",       "#1F78B4",
        "sim_hi_r", "Simulation Contact high", "#A6CEE3",
        "fb_lo_r",  "Data Contact low",        "#FF7F00",
        "sim_lo_r", "Simulation Contact low",  "#FDBF6F"
    )

    slw <- if ("mean" %in% plot_opts) 0.7 else 0.2
    slt <- if ("mean" %in% plot_opts) "dashed" else "solid"

    description <- params$description |>
        str_split_1(", ") |>
        str_subset("convergence|coverage|pedigree|GRM", negate = TRUE) |>
        str_flatten_comma()


    l2p <- 1 / ggplot2::.pt

    sim_vals <- if ("ribbon" %notin% plot_opts) {
        geom_ribbon(aes(x = time, ymin = hdi1, ymax = hdi2,
                group = gp, fill = str),
            data_sim,
            alpha = 0.5,
            show.legend = FALSE)
    } else {
        geom_line(aes(x = time, y = survival,
                      group = gp, colour = str),
                  data_svm[src == "sim"],
                  linewidth = 0.5 * l2p)
    }

    sim_means <- if (all(c("ribbon", "mean") %in% plot_opts)) {
        geom_line(aes(x = time, y = survival,
                      group = gp, colour = str),
                  linewidth = 1 * l2p,
                  linetype = "dashed",
                  data_sim,
                  alpha = 0.5,
                  show.legend = FALSE)
    }

    fb_vals <- geom_line(aes(x = time, y = survival,
                             group = gp, colour = str),
                         data_svm[src == "fb"])

    plt <- ggplot() +
        sim_vals +
        sim_means +
        fb_vals +
        scale_fill_manual(breaks = scm$breaks,
                          labels = scm$labels,
                          values = scm$values) +
        scale_colour_manual(breaks = scm$breaks,
                            labels = scm$labels,
                            values = scm$values) +
        lims(y = c(0, 100)) +
        labs(colour = "Source",
             x = "Time (days)",
             y = "Survival (%)",
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

    list(plt = plt, data = data_svm)
}
