library(stringr)
library(ggplot2)
library(viridisLite)

plot_model <- function(popn, params) {
    pm <- get(str_c("plot_", str_remove(params$model_type, "_res")))

    pm(popn, params)
}


plot_SxxDR <- function(popn, params) {
    message("Plotting SxxDR model")

    N <- popn[sdp == "progeny", .N]
    tmax <- params$tmax

    events <- make_time_series_seidr(popn, params)
    events[, S := S + E + I]
    events[, c("E", "I", "ID") := NULL]
    # e2 <- melt(events, id = "time")

    if (params$sim_new_data == "no") {
        events <- events[time != max(time)]
    }
    
    plt <- ggplot(events, aes(x = time)) +
        geom_line(aes(y = S / N, colour = "Susceptible?"), linewidth = 1.2) +
        geom_line(aes(y = D / N, colour = "Detectable"),   linewidth = 0.6) +
        geom_line(aes(y = R / N, colour = "Removed"),      linewidth = 1.2) +
        scale_colour_manual("Compartments",
                            breaks = c("Susceptible?", "Detectable", "Removed"),
                            values = c("blue", "red", "green")) +
        labs(title = "SxxDR", x = "Time (days)", y = "Proportion") +
        coord_cartesian(xlim = c(0, min(tmax, max(events$time), na.rm = TRUE))) +
        theme_bw() +
        theme(legend.position = "bottom")

    plt
}

plot_SEIDR <- function(popn, params) {
    message("Plotting SEIDR model")

    N <- popn[sdp == "progeny", .N]
    tmax <- max(params$tmax)

    events <- make_time_series_seidr(popn, params)
    # e2 <- melt(events, id = "time")

    plt <- ggplot(events, aes(x = time)) +
        geom_line(aes(y = S / N,  colour = "Susceptible"),  linewidth = 1.2) +
        geom_line(aes(y = E / N,  colour = "Exposed"),      linewidth = 1.2) +
        geom_line(aes(y = I / N,  colour = "Undetectable"), linewidth = 0.6) +
        geom_line(aes(y = D / N,  colour = "Detectable"),   linewidth = 0.6) +
        geom_line(aes(y = ID / N, colour = "Infectious"),   linewidth = 1.2) +
        geom_line(aes(y = R / N,  colour = "Removed"),      linewidth = 1.2) +
        scale_colour_manual("Compartments",
                            breaks = c("Susceptible", "Exposed", "Undetectable",
                                       "Detectable", "Infectious", "Removed"),
                            values = c(viridis(5), "red")) +
        labs(title = "SEIDR", x = "Time (days)", y = "Proportion") +
        coord_cartesian(xlim = c(0, min(tmax, max(events$time), na.rm = TRUE))) +
        theme_bw() +
        theme(legend.position = "bottom")

    if (params$show_plots) {
        print(plt)
    }
    plt
}


plot_SIDR <- function(popn, params) {
    message("Plotting SIDR model")

    N <- popn[sdp == "progeny", .N]
    tmax <- params$tmax

    events <- make_time_series_sidr(popn, params)
    # e2 <- melt(events, id = "time")

    plt <- ggplot(events, aes(x = time)) +
        geom_line(aes(y = S / N,  colour = "Susceptible"),  linewidth = 1.2) +
        geom_line(aes(y = I / N,  colour = "Undetectable"), linewidth = 0.6) +
        geom_line(aes(y = D / N,  colour = "Detectable"),   linewidth = 0.6) +
        geom_line(aes(y = ID / N, colour = "Infectious"),   linewidth = 1.2) +
        geom_line(aes(y = R / N,  colour = "Removed"),      linewidth = 1.2) +
        scale_colour_manual("Compartments",
                            breaks = c("Susceptible", "Undetectable", "Detectable", "Infectious", "Removed"),
                            values = c("blue", "mediumpurple", "purple", "red", "green")) +
        labs(title = "SIDR", x = "Time (days)", y = "Proportion") +
        coord_cartesian(xlim = c(0, min(tmax, max(events$time), na.rm = TRUE))) +
        theme_bw() +
        theme(legend.position = "bottom")

    print(plt)
}


plot_SEIR <- function(popn, params) {
    message("Plotting SEIR model")

    N <- popn[sdp == "progeny", .N]
    tmax <- params$tmax

    events <- make_time_series_seir(popn, params)
    # e2 <- melt(events, id = "time")

    plt <- ggplot(events, aes(x = time)) +
        geom_line(aes(y = S / N, colour = "Susceptible"), linewidth = 1.2) +
        geom_line(aes(y = E / N, colour = "Exposed"),     linewidth = 1.2) +
        geom_line(aes(y = I / N, colour = "Infectious"),  linewidth = 1.2) +
        geom_line(aes(y = R / N, colour = "Removed"),     linewidth = 1.2) +
        scale_colour_manual("Compartments",
                            breaks = c("Susceptible", "Exposed", "Infectious", "Removed"),
                            values = c("blue", "pink", "red", "green")) +
        labs(title = "SEIR", x = "Time (days)", y = "Proportion") +
        coord_cartesian(xlim = c(0, min(tmax, max(events$time), na.rm = TRUE))) +
        theme_bw() +
        theme(legend.position = "bottom")

    print(plt)
}


plot_SIR <- function(popn, params) {
    message("Plotting SIR model")

    N <- popn[sdp == "progeny", .N]
    tmax <- params$tmax

    events <- make_time_series_sir(popn, params)

    plt <- ggplot(events, aes(x = time)) +
        # geom_line(linewidth = 1.2) +
        geom_line(aes(y = S / N, colour = "Susceptible"), linewidth = 1.2) +
        geom_line(aes(y = I / N, colour = "Infectious"),  linewidth = 1.2) +
        geom_line(aes(y = R / N, colour = "Removed"),     linewidth = 1.2) +
        scale_colour_manual("Compartments",
                            breaks = c("Susceptible", "Infectious", "Removed"),
                            values = c("blue", "red", "green")) +
        labs(title = "SIR", x = "Time (days)", y = "Proportion") +
        coord_cartesian(xlim = c(0, min(tmax, max(events$time), na.rm = TRUE))) +
        theme_bw() +
        theme(legend.position = "bottom")

    print(plt)
}


plot_SIS <- function(popn, params) {
    message("Plotting SIS model - warning: this is currently broken")

    N <- popn[sdp == "progeny", .N]
    tmax <- params$tmax

    events <- make_time_series_sis(popn, params)

    plt <- ggplot(events, aes(x = time)) +
        # geom_line(linewidth = 1.2) +
        geom_line(aes(y = S / N, colour = "Susceptible"), linewidth = 1.2) +
        geom_line(aes(y = I / N, colour = "Infectious"),  linewidth = 1.2) +
        scale_colour_manual("Compartments",
                            breaks = c("Susceptible", "Infectious"),
                            values = c("blue", "red")) +
        labs(title = "SIR", x = "Time (days)", y = "Proportion") +
        coord_cartesian(xlim = c(0, min(tmax, max(events$time), na.rm = TRUE))) +
        theme_bw() +
        theme(legend.position = "bottom")

    print(plt)
}
