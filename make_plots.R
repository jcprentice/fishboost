library(ggplot2)

plot_model <- function(pop, params) {
    pm <- get(paste0("plot_", sub("_res", "", params$model_type)))

    pm(pop, params)
}


plot_SxxDR <- function(pop, params) {
    message("Plotting SxxDR model")

    N <- pop[sdp == "progeny", .N]
    tmax <- params$tmax

    events <- make_time_series_seidr(pop, params)
    events[, S := S + E + I]
    events[, c("E", "I", "ID") := NULL]
    # e2 <- melt(events, id = "time")

    if (params$use_fb_data) {
        events <- events[time != max(time)]
    }

    plt <- ggplot(events, aes(x = time)) +
        geom_line(aes(y = S / N, colour = "Susceptible?"), linewidth = 1.2) +
        geom_line(aes(y = D / N, colour = "Detectable"),   linewidth = 0.6) +
        geom_line(aes(y = R / N, colour = "Recovered"),    linewidth = 1.2) +
        scale_colour_manual("Compartments",
                            breaks = c("Susceptible?", "Detectable", "Recovered"),
                            values = c("blue", "red", "green")) +
        labs(title = "SxxDR", x = "Time (days)", y = "Proportion") +
        coord_cartesian(xlim = c(0, min(tmax, max(events$time), na.rm = TRUE))) +
        theme(legend.position = "bottom")

    plt
}

plot_SEIDR <- function(pop, params) {
    message("Plotting SEIDR model")

    N <- pop[sdp == "progeny", .N]
    tmax <- params$tmax

    events <- make_time_series_seidr(pop, params)
    # e2 <- melt(events, id = "time")

    plt <- ggplot(events, aes(x = time)) +
        geom_line(aes(y = S / N, colour = "Susceptible"),  linewidth = 1.2) +
        geom_line(aes(y = E / N, colour = "Exposed"),      linewidth = 1.2) +
        geom_line(aes(y = I / N, colour = "Undetectable"), linewidth = 0.6) +
        geom_line(aes(y = D / N, colour = "Detectable"),   linewidth = 0.6) +
        geom_line(aes(y = ID / N, colour = "Infectious"),  linewidth = 1.2) +
        geom_line(aes(y = R / N, colour = "Recovered"),    linewidth = 1.2) +
        scale_colour_manual("Compartments",
                            breaks = c("Susceptible", "Exposed", "Undetectable", "Detectable", "Infectious", "Recovered"),
                            values = c("blue", "pink", "mediumpurple", "purple", "red", "green")) +
        labs(title = "SEIDR", x = "Time (days)", y = "Proportion") +
        coord_cartesian(xlim = c(0, min(tmax, max(events$time), na.rm = TRUE))) +
        theme(legend.position = "bottom")

    if (params$show_plots) {
        print(plt)
    }
    plt
}


plot_SIDR <- function(pop, params) {
    message("Plotting SIDR model")

    N <- pop[sdp == "progeny", .N]
    tmax <- params$tmax

    events <- make_time_series_sidr(pop, params)
    # e2 <- melt(events, id = "time")

    plt <- ggplot(events, aes(x = time)) +
        geom_line(aes(y = S / N, colour = "Susceptible"),  linewidth = 1.2) +
        geom_line(aes(y = I / N, colour = "Undetectable"), linewidth = 0.6) +
        geom_line(aes(y = D / N, colour = "Detectable"),   linewidth = 0.6) +
        geom_line(aes(y = ID / N, colour = "Infectious"),  linewidth = 1.2) +
        geom_line(aes(y = R / N, colour = "Recovered"),    linewidth = 1.2) +
        scale_colour_manual("Compartments",
                            breaks = c("Susceptible", "Undetectable", "Detectable", "Infectious", "Recovered"),
                            values = c("blue", "mediumpurple", "purple", "red", "green")) +
        labs(title = "SIDR", x = "Time (days)", y = "Proportion") +
        coord_cartesian(xlim = c(0, min(tmax, max(events$time), na.rm = TRUE))) +
        theme(legend.position = "bottom")

    print(plt)
}


plot_SEIR <- function(pop, params) {
    message("Plotting SEIR model")

    N <- pop[sdp == "progeny", .N]
    tmax <- params$tmax

    events <- make_time_series_seir(pop, params)
    # e2 <- melt(events, id = "time")

    plt <- ggplot(events, aes(x = time)) +
        geom_line(aes(y = S / N, colour = "Susceptible"), linewidth = 1.2) +
        geom_line(aes(y = E / N, colour = "Exposed"),     linewidth = 1.2) +
        geom_line(aes(y = I / N, colour = "Infectious"),  linewidth = 1.2) +
        geom_line(aes(y = R / N, colour = "Recovered"),   linewidth = 1.2) +
        scale_colour_manual("Compartments",
                            breaks = c("Susceptible", "Exposed", "Infectious", "Recovered"),
                            values = c("blue", "pink", "red", "green")) +
        labs(title = "SEIR", x = "Time (days)", y = "Proportion") +
        coord_cartesian(xlim = c(0, min(tmax, max(events$time), na.rm = TRUE))) +
        theme(legend.position = "bottom")

    print(plt)
}


plot_SIR <- function(pop, params) {
    message("Plotting SIR model")

    N <- pop[sdp == "progeny", .N]
    tmax <- params$tmax

    events <- make_time_series_sir(pop, params)

    plt <- ggplot(events, aes(x = time)) +
        # geom_line(linewidth = 1.2) +
        geom_line(aes(y = S / N, colour = "Susceptible"), linewidth = 1.2) +
        geom_line(aes(y = I / N, colour = "Infectious"),  linewidth = 1.2) +
        geom_line(aes(y = R / N, colour = "Recovered"),   linewidth = 1.2) +
        scale_colour_manual("Compartments",
                            breaks = c("Susceptible", "Infectious", "Recovered"),
                            values = c("blue", "red", "green")) +
        labs(title = "SIR", x = "Time (days)", y = "Proportion") +
        coord_cartesian(xlim = c(0, min(tmax, max(events$time), na.rm = TRUE))) +
        theme(legend.position = "bottom")

    print(plt)
}


plot_SIS <- function(pop, params) {
    message("Plotting SIS model - warning: this is currently broken")

    N <- pop[sdp == "progeny", .N]
    tmax <- params$tmax

    events <- make_time_series_sis(pop, params)

    plt <- ggplot(events, aes(x = time)) +
        # geom_line(linewidth = 1.2) +
        geom_line(aes(y = S / N, colour = "Susceptible"), linewidth = 1.2) +
        geom_line(aes(y = I / N, colour = "Infectious"),  linewidth = 1.2) +
        scale_colour_manual("Compartments",
                            breaks = c("Susceptible", "Infectious"),
                            values = c("blue", "red")) +
        labs(title = "SIR", x = "Time (days)", y = "Proportion") +
        coord_cartesian(xlim = c(0, min(tmax, max(events$time), na.rm = TRUE))) +
        theme(legend.position = "bottom")

    print(plt)
}
