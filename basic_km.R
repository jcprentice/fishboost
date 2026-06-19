{
    library(data.table)
    library(stringr)
    library(ggplot2)
}

#' Generate a basic KM plot from a single data.table
#'
#' @param popn A population data.table
#' @param params A parameters list
#'
#' @returns A KM plot

basic_km <- function(popn, params) {
    show_plots <- params$show_plots

    if (FALSE) {
        popn <- readRDS("fb_data/fb_12_rpw.rds")
        show_plots <- TRUE

        popn[Tdeath %in% c(104, 160), Tdeath := NA]
    }


    cols <- intersect(c("sire", "dam", "trial", "donor", "Tinf", "Tsign", "Tdeath"),
                      names(popn))

    x <- popn[sdp == "progeny", ..cols]

    x[, family := .GRP, .(sire, dam, trial)] |>
        setcolorder("family")

    # Drop small families
    x[, N := .N, family]
    x <- x[N >= 10]
    x[, N := NULL]

    if ("Tsign" %notin% names(x)) {
        x[, Tsign := Tinf]
    }

    # If Tsign is missing but Tdeath is not, then let Tsign = Tdeath, otherwise
    # set to tmax + 1 (lines should extend past the edge of the figure)
    x[is.na(Tsign), Tsign := fifelse(!is.na(Tdeath), Tdeath, tmax[trial] + 1)]
    # RP can be extended from Tsign to tmax + 1 if Tdeath is missing
    x[, RP := Tdeath - Tsign]
    x[is.na(RP), RP := tmax[trial] - Tsign]
    x[, RP := pmax(RP, 0)]

    sv_curve <- function(x) c(0, sort(x, na.last = TRUE))

    x1 <- x[, .(donor = fifelse(any(donor == 1), "donor", "recip"),
                survival = seq(100, 0, length.out = .N + 1),
                Tsign = sv_curve(Tsign),
                RP    = sv_curve(RP)),
            .(family, trial)] |>
        melt(measure.vars = c("Tsign", "RP"),
             value.name = "time")

    # Resample at t = 0,1,...,tmax
    x2 <- x1[, {
        times <- seq(0, tmax[trial[[1]]])
        list(time = times,
             survival = approx(time, survival, times,
                               method = "constant",
                               ties = list("ordered", max))$y)
    }, .(family, trial, donor, variable)]

    plt <- ggplot(x2,
                  aes(x = time, y = survival,
                      colour = donor, group = family)) +
        geom_line() +
        # geom_point(size = 0.5) +
        scale_colour_manual("Inoculated Family",
                            breaks = c("donor", "recip"),
                            labels = c("Yes", "No"),
                            values = c("#33A02C", "#1F78B4")) +
        lims(y = c(0, 100)) +
        labs(x = "Time (days)",
             y = "Survival (%)") +
        facet_grid(rows = vars(trial),
                   cols = vars(variable),
                   scales = "free_x",
                   labeller = labeller(
                       variable = c(Tinf  = "Time of infection",
                                    Tsign = "Time to first visual signs",
                                    RP    = "Time from visual signs to death"),
                       trial = c("1" = "Trial 1",
                                 "2" = "Trial 2"))) +
        theme_bw() +
        theme(text = element_text(size = 5),
              plot.title = element_markdown(size = 7),
              strip.background = element_blank(),
              # panel.background = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              legend.position = "bottom")
    plt
}
