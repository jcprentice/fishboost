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


    cols <- intersect(c("sire", "trial", "donor", "Tinf", "Tsym", "Tdeath"),
                      names(popn))

    x <- popn[sdp == "progeny", ..cols]

    x[, N := .N, .(trial, sire)]
    x <- x[N >= 10]
    x[, N := NULL]

    if ("Tsym" %notin% names(x)) {
        x[, Tsym := Tinf]
    }

    x[, RP := Tdeath - Tsym]
    x[is.na(Tsym), Tsym := Tdeath]


    x1 <- x[, .(donor = fifelse(any(donor == 1), "donor", "recip"),
                Tsym = c(0, sort(Tsym, na.last = TRUE)),
                RP   = c(0, sort(RP,   na.last = TRUE))),
            .(sire, trial)]
    x1[, grp := .GRP, .(sire, trial)]
    x1[, survival := seq(1, 0, length.out = .N), grp]

    x2 <- melt(x1, measure.vars = c("Tsym", "RP"),
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
                       variable = c(Tinf = "Proportion of family uninfected vs time",
                                    Tsym = "Proportion of family with no symptoms vs time",
                                    RP   = "Proportion of family surviving vs time"),
                       trial = c("1" = "Trial 1",
                                 "2" = "Trial 2"))) +
        theme_bw() +
        theme(panel.background = element_blank(),
              # panel.grid.major = element_blank(),
              # panel.grid.minor = element_blank(),
              legend.position = "bottom")
    plt
}
