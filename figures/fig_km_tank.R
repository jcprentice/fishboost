{
    library(data.table)
    library(stringr)
    library(purrr)
    library(HDInterval)
    library(ggplot2)
    library(ggtext)
    library(cowplot)

    source("figures/theme_natcom.R")
}

fig_km_tank <- function(dataset = "fb-test",
                        scen = 2,
                        plot_opts = c(),
                        title_str = NULL,
                        DEBUG = FALSE) {

    if (FALSE) {
        dataset <- "fb-test"
        scen <- 7
        plot_opts <- c("fb_only", "extremes", "t1", "t2", "t2_width", "mean", "ribbon")
        plot_opts <- c("Tsign", "t1", "ribbon")
        title_str <- "Tsign<br>Multi-stage endurance"
        DEBUG <- TRUE
    }

    {
        km_data  <- readRDS(str_glue("datasets/{dataset}/meta/km_data_ps.rds"))
        data     <- copy(km_data[[scen]]$data)
        params   <- km_data[[scen]]$params
        opts     <- km_data[[scen]]$opts
    }

    if ("fb_only" %in% plot_opts) {
        if (DEBUG) message("- Dropping all simulated values")
        data <- data[src == "fb"]
    }

    # One and only one of these should be present in plot_opts
    timings <- c("Tinf", "Tsign", "RP")

    if (length(intersect(timings, plot_opts)) != 1) {
        message("plot_opts: ", str_flatten_comma(plot_opts))
        stop("Only one timing should be present")
    }

    # Filter by trial
    if ("t1" %in% plot_opts & "t2" %notin% plot_opts) {
        if (DEBUG) message(" - Keeping Trial 1 only")
        data <- data[trial == 1]
    } else if ("t1" %notin% plot_opts & "t2" %in% plot_opts) {
        if (DEBUG) message(" - Keeping Trial 2 only")
        data <- data[trial == 2]
    }


    # If Tsign is missing but Tdeath is not, then let Tsign = Tdeath, otherwise
    # set to tmax + 10 (lines should extend past the edge of the figure).
    # RP can be extended from Tsign to tmax + 1 if Tdeath is missing.
    {
        tmax <- params$tmax
        data[is.na(Tsign), Tsign := fifelse(!is.na(Tdeath), Tdeath, tmax[trial] + 1)]
        data[, RP := Tdeath - Tsign]
        data[is.na(RP), RP := tmax[trial] - Tsign]
        data[, RP := pmax(RP, 0)]
    }

    # Make a copy of the data and add in rows representing t=0, sort times, make
    # sure to include NA values last
    sv_curve <- function(x) c(0, sort(x, na.last = TRUE))

    data_sv <- data[, .(survival = seq(100, 0, length.out = .N + 1),
                        Tinf  = sv_curve(Tinf),
                        Tsign = sv_curve(Tsign),
                        RP    = sv_curve(RP)),
                    .(id, src, group, trial)] |>
        setorder(id, group, trial)

    # All lines in the plot will be based on this grouping
    data_sv[, gp := .GRP, .(id, group, trial)]

    # Set extreme families
    if ("extremes" %in% plot_opts) {
        if (DEBUG) message("- Using extreme groups")
        foo <- data_sv[src == "fb", .(group, trial, Tsign)]
        # Filter out any sires with fewer than 5 non-NA values
        foo[, p := sum(!is.na(Tsign)), .(group, trial)]
        foo1 <- foo[p > 5,
                    .(mu = mean(Tsign, na.rm = TRUE)),
                    .(group, trial)] |>
            setorder(group, trial, mu)
        ids <- foo1[, .(group = group[c(which.min(mu), which.max(mu))],
                        es = c("lo", "hi")),
                    trial]
        data_sv <- merge(data_sv, ids, by = c("trial", "group"))
    } else {
        if (DEBUG) message("- Using all families")
        data_sv[, es := "X"]
    }

    # Melt so we can use facet_wrap, and resample times
    tmp <- melt(data_sv,
                measure.vars = patterns("Tinf|Tsign|RP"),
                value.name = "time")

    # Use approxfun to convert to t = 0, 1, ..., tmax

    data_svm <- tmp |>
        _[variable %in% plot_opts, {
            times <- seq(0, tmax[trial])
            list(time = times,
                 survival = approx(time, survival, times,
                                   method = "constant",
                                   ties = list("ordered", max))$y)
        },
        .(id, group, trial, src, es, gp, variable)]

    data_svm[, `:=`(survival = nafill(survival, type = "locf"),
                    str = str_c(src, es, sep = "_") |> str_remove("_X"))] |>
        setcolorder("str", after = "es")


    # Get simulation means and HDI
    data_sim <- data_svm[src == "sim",
                         .(gp = first(gp),
                           survival = mean(survival),
                           hdi1 = hdi(survival)[["lower"]],
                           hdi2 = hdi(survival)[["upper"]]),
                         .(trial, src, es, str, variable, time)] |>
        setcolorder(names(data_svm), skip_absent = TRUE)


    scm <- rowwiseDT(
        breaks=,    labels=,                      values=,
        "fb",     "Data Contact",               "#1F78B4",
        "sim",    "Simulation Contact",         "#A6CEE3",
        "fb_hi",  "Data Contact highest",       "#1F78B4",
        "sim_hi", "Simulation Contact highest", "#A6CEE3",
        "fb_lo",  "Data Contact lowest",        "#FF7F00",
        "sim_lo", "Simulation Contact lowest",  "#FDBF6F"
    )

    description <- params$description |>
        str_split_1(", ") |>
        str_subset("convergence|coverage|pedigree|GRM", negate = TRUE) |>
        str_flatten_comma()


    # Plot ----

    l2p <- 1 / ggplot2::.pt

    sim_vals <- if ("ribbon" %in% plot_opts) {
        geom_ribbon(aes(x = time, ymin = hdi1, ymax = hdi2,
                        group = gp, fill = str),
                    data_sim,
                    linewidth = 0.5 * l2p,
                    alpha = 0.5)
    } else {
        geom_line(aes(x = time, y = survival,
                      group = gp, colour = str),
                  data_svm[src == "sim"])
    }

    sim_means <- if ("mean" %in% plot_opts) {
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

    Tmax <- if ("t2_width" %in% plot_opts) 160 else max(data_svm$time, na.rm = TRUE)

    plt <- ggplot() +
        sim_vals +
        sim_means +
        fb_vals +
        scale_fill_manual(breaks = scm$breaks,
                          labels = scm$labels,
                          values = scm$values,
                          drop = FALSE) +
        scale_colour_manual(breaks = scm$breaks,
                            labels = scm$labels,
                            values = scm$values,
                            drop = FALSE) +
        # guides(fill   = guide_legend(ncol = 1),
        #        colour = guide_legend(ncol = 1)) +
        lims(x = c(0, Tmax),
             y = c(0, 100)) +
        labs(fill = NULL,
             colour = NULL,
             x = "Time (days)",
             y = "Survival (%)",
             title = title_str) +
        theme_natcom()
    plt
}

{
    p1 <- fig_km_tank("fb-test", 2, c("Tsign", "ribbon"),
                        "Time to visual signs<br>Multi-stage endurance")

    p2 <- fig_km_tank("fb-test", 2, c("Tsign", "extremes", "ribbon"),
                        "Time to visual signs, extremes<br>Multi-stage endurance")

    p3 <- fig_km_tank("fb-test", 3, c("Tsign", "extremes", "ribbon"),
                        "Time to visual signs, extremes<br>No genetic variance")

    p4 <- fig_km_tank("fb-test", 2, c("RP", "ribbon"),
                        "Time from visual signs to death<br>Multi-stage endurance")

    p5 <- fig_km_tank("fb-test", 2, c("RP", "extremes", "ribbon"),
                        "Time from visual signs to death, extremes<br>Multi-stage endurance")

    p6 <- fig_km_tank("fb-test", 3, c("RP", "extremes", "ribbon"),
                        "Time from visual signs to death, extremes<br>No genetic variance")

    scm <- rowwiseDT(
        breaks=,    labels=,                      values=,
        "fb",     "Data Contact",               "#1F78B4",
        "sim",    "Simulation Contact",         "#A6CEE3",
        "fb_hi",  "Data Contact highest",       "#1F78B4",
        "sim_hi", "Simulation Contact highest", "#A6CEE3",
        "fb_lo",  "Data Contact lowest",        "#FF7F00",
        "sim_lo", "Simulation Contact lowest",  "#FDBF6F"
    )

    N <- nrow(scm)
    p_leg <- ggplot(data.table(x = rep(seq_len(N), 10),
                               y = runif(10 * N),
                               z = rep(scm$breaks, each = 10))) +
        aes(x, y, group = z, colour = z) +
        geom_line() +
        guides(colour = guide_legend(nrow = 2)) +
        labs(colour = NULL) +
        scale_colour_manual(breaks = scm$breaks,
                            labels = scm$labels,
                            values = scm$values) +
        theme_classic() +
        theme(text = element_text(size = 5))
    p_leg
    p_legend <- get_legend(p_leg)
}

plt1 <- plot_grid(p1 + theme(legend.position = "none"),
                  p2 + theme(legend.position = "none"),
                  p3 + theme(legend.position = "none"),
                  p4 + theme(legend.position = "none"),
                  p5 + theme(legend.position = "none"),
                  p6 + theme(legend.position = "none"),
                  labels = "auto",
                  label_size = 10,
                  align = "hv", ncol = 3)

plt <- plot_grid(plt1, p_legend,
                 ncol = 1, rel_heights = c(1, 0.08))
plt

ggsave("gfx/km_plots-tank.pdf", plt,
       width = 18.3, height = 17, units = "cm", dpi = "print")
ggsave("gfx/km_plots-tank.png", plt,
       width = 18.3, height = 17, units = "cm", dpi = "print")
