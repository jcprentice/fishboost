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
    cols <- intersect(c("sire", "trial", "Tinf", "Tsym", "Tdeath"),
                      names(popn))
    
    x <- popn[sdp == "progeny", ..cols]
    
    if ("Tsym" %notin% names(x)) {
        x[, Tsym := Tinf]
    }
    
    x[, RP := Tdeath - Tsym]
    
    x1 <- x[, .(Tsym = c(0, sort(Tsym, na.last = TRUE)),
                RP   = c(0, sort(RP,   na.last = TRUE))),
            .(sire, trial)]
    x1[, grp := .GRP, .(sire, trial)]
    x1[, survival := seq(1, 0, length.out = .N), grp]
    
    x2 <- melt(x1, measure.vars = c("Tsym", "RP"),
               value.name = "time")
    
    plt <- ggplot(x2, aes(x = time, y = survival, group = grp)) +
        geom_line(colour = "red") +
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
        theme_bw()
    
    if (params$show_plots) print(plt)
    
    plt
}
