{
    library(data.table)
    library(stringr)
    library(purrr)
    library(ggplot2)
}

plot_km_variant <- function(data_list,
                            KM_measure = "Tsym",
                            KM_var = "sire",
                            trials = 12) {
    
    # data_list <- km_data[[1]]; KM_measure <- "Tsym"; KM_var <- "sire"; trials <- 12

    # Extract parts from data_list
    params       <- data_list$params
    dataset      <- data_list$opts$dataset
    scenario     <- data_list$opts$scenario
    apply_gammas <- data_list$opts$apply_gammas
    n_plots      <- data_list$opts$n_plots

    # Make sure to copy this otherwise we'll end up modifying it
    data         <- copy(data_list$data)

    # Rename the 2 parts we want to plot
    setnames(data, c(KM_var, KM_measure), c("variable", "KM"))

    # Make a copy of the data and add in rows representing t=0. Put NA's last
    # and set them to Inf
    data_t0 <- data[, .(KM = c(0, sort(KM, na.last = TRUE)),
                        src = src[c(1, 1:.N)]),
                    .(id, variable, trial)]
    # data_t0[is.na(KM), KM := Inf]

    # Rename src
    data_t0[, src := str_c(src, trial)]

    # Split by trial
    if (trials %in% 1:2) {
        data_t0 <- data_t0[trial == trials]
    }


    # Survival curves by id and sire
    data_t0[, survival := seq(1, 0, length.out = .N), .(id, variable)]

    # Cut off the last 10% since it tends to dominate... or just stop at 160 days
    # Tmax <- data_t0[!is.na(Tdeath), quantile(Tdeath, 0.9)[[1]]]
    Tmax <- data_t0[id == max(id), max(KM, na.rm = TRUE)]

    # Create a column to group by
    data_t0[, gp := factor(str_c(sprintf("%02d", id), "_", sprintf("%02d", variable), "_", trial))]

    # String interpolation for the title
    wt_str <- if (apply_gammas) "gamma" else "exponential"

    KMm_str <- switch(KM_measure,
                      "Tsym" = "Time to symptoms",
                      "Tdeath" = "Time to death",
                      "RP" = "Survival after symptoms",
                      "Measure")
    KMv_str <- switch(KM_var,
                      "sire" = "family",
                      "group" = "tank",
                      "group")
    KM_xlim <- switch(KM_measure,
                      "Tsym" = 300,
                      "Tdeath" = 300,
                      "RP" = 150,
                      max(data_t0$KM))


    # Plot
    plt <- ggplot(data_t0) +
        # geom_step
        geom_line(aes(x = KM, y = survival, group = gp, colour = src)) +
        # coord_cartesian(xlim = c(0, min(200, max(data_t0$KM)))) +
        lims(x = c(0, 200), y = 0:1) +
        labs(x = str_glue("{KMm_str} (days)"),
             y = "Survival",
             title = str_glue("{KMm_str} by {KMv_str}"),
             # title = str_c("KM plot for FB vs simulated data, parasites, (", wt_str, " waiting times)"),
             colour = "Source") +
        scale_colour_manual(breaks = c("fb1", "fb2", "sim1", "sim2"),
                            labels = c("FB Trial 1", "FB Trial 2",
                                       "Sim Trial 1", "Sim Trial 2"),
                            values = c("blue", "red", "lightblue", "pink")) +
        theme(panel.background = element_blank())
    plt



    list(plt = plt, data = data_t0)
}


