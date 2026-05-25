library(data.table)
library(stringr)
library(ggplot2)

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


    cols <- intersect(c("sire", "trial", "donor", "Tinf", "Tsign", "Tdeath"),
                      names(popn))

    x <- popn[sdp == "progeny", ..cols]

    x[, N := .N, .(trial, sire)]
    x <- x[N >= 10]
    x[, N := NULL]

    if ("Tsign" %notin% names(x)) {
        x[, Tsign := Tinf]
    }

    x[, RP := Tdeath - Tsign]
    x[is.na(Tsign), Tsign := Tdeath]

    sv_curve <- function(x) c(0, sort(x, na.last = TRUE))

    x1 <- x[, .(donor = fifelse(any(donor == 1), "donor", "recip"),
                Tsign = sv_curve(Tsign),
                RP    = sv_curve(RP)),
            .(sire, trial)]
    x1[, grp := .GRP, .(sire, trial)]
    x1[, survival := seq(1, 0, length.out = .N), grp]

    x2 <- melt(x1, measure.vars = c("Tsign", "RP"),
               value.name = "time")

    plt <- ggplot(x2, aes(x = time, y = survival, colour = donor, group = grp)) +
        geom_line() +
        geom_point(size = 0.5) +
        scale_colour_manual("Inoculated",
                            breaks = c("donor", "recip"),
                            labels = c("Yes", "No"),
                            values = c("#FF7F00", "#1F78B4")) +
        labs(x = "Time (days)",
             y = "Survival") +
        facet_grid(rows = vars(trial),
                   cols = vars(variable),
                   scales = "free_x",
                   labeller = labeller(
                       variable = c(Tinf  = "Time of infection",
                                    Tsign = "Time to first signs",
                                    RP    = "Period from signs to death"),
                       trial = c("1" = "Trial 1",
                                 "2" = "Trial 2"))) +
        theme_bw() +
        theme(panel.background = element_blank(),
              # panel.grid.major = element_blank(),
              # panel.grid.minor = element_blank(),
              legend.position = "bottom")
    plt
}
