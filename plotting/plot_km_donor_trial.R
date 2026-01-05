{
    library(data.table)
    library(stringr)
    library(purrr)
    library(ggplot2)
}

plot_km_donor_trial <- function(data_list, trials = 12, plotopts = c()) {
    
    if (FALSE) {
        i <- 2
        data <- copy(km_data[[i]]$data)
        params <- km_data[[i]]$params
        opts <- km_data[[i]]$opts
        trials <- 12
        plotopts = c("drop_small_groups")
    } else {
        data   <- copy(data_list$data)
        params <- data_list$params
        opts   <- data_list$opts
    }    
    
    description <- params$description |>
        str_remove_all(", convergence|, coverage|, pedigree|, GRM \\w*") |>
        str_replace_all(c("inf_model 1" = "inf: I = D",
                          "inf_model 2" = "inf: I = 0.1*D",
                          "inf_model 3" = "inf: Don = 0.1*Rec",
                          "inf_model 4" = "inf: Don = r*Rec")) 
    
    # Choose between actual Tinfs and BICI's inferred values for FB data
    if ("use_sire_Tinfs" %in% plotopts && "Tinf_sire" %in% names(data)) {
        fb_Tinfs <- with(params, get_median_Tinfs(dataset, scenario, replicate))
        
        N1 <- data[src == "fb", .N]
        N2 <- nrow(fb_Tinfs)
        data[src == "fb", Tinf2 := fb_Tinfs$Tinf[seq(1 + N2 - N1, N2)]]
    }
    
    # If Tsym is missing but Tdeath is not, then let Tsym = Tdeath
    data[, Tsym2 := fifelse(is.na(Tsym) & !is.na(Tdeath), Tdeath, Tsym)]
    
    # Make a copy of the data and add in rows representing t=0, sort times
    data_t0 <- data[, .(Tinf = c(0, sort(Tinf,  na.last = TRUE)),
                        Tsym = c(0, sort(Tsym2, na.last = TRUE)),
                        RP   = c(0, sort(RP,    na.last = TRUE)),
                        src = src[c(1, 1:.N)]),
                    .(id, donor, trial)]
    
    setorder(data_t0, id, donor, trial)
    
    # Rename src
    data_t0[, src := str_c(src, "_", fifelse(donor == 1, "d", "r"))]
    
    # Create a column to group by, ensure this is wide enough to handle all ids
    # data_t0[, gp := factor(str_c(format(id), "_", donor, "_", trial))]
    data_t0[, gp := .GRP, .(id, donor, trial)]
    
    # Split by trial
    if (trials %in% 1:2) {
        data_t0 <- data_t0[trial == trials]
    }
    
    # Survival curves by id and sire
    data_t0[, survival := seq(1, 0, length.out = .N), gp]
    
    # Melt so we can use facet_wrap
    data_t1 <- melt(data_t0,
                    measure.vars = c("Tinf", "Tsym", "RP"))
    
    if ("show_Tinfs" %notin% plotopts) {
        data_t1 <- data_t1[variable != "Tinf"]
    }
    
    # Drop the really high values
    # data_t1[variable == "Tsym" & value > 200, value := NA]
    # data_t1[variable == "RP"   & value > 160, value := NA]
    # data_t1[value > 160, value := NA]
    # data_t1[(trial == 1 & value >= 104) | (trial == 2 & value >= 160), value := NA]
    maxvals <- c(104, 160)
    data_t1[, value := fifelse(value >= maxvals[trial], maxvals[trial], value)]
    
    um <- if (opts$use_means) ", mean" else ", samples"
    
    # Plot
    plt <- ggplot(data_t1) +
        # geom_step
        geom_line(aes(x = value,
                      y = survival,
                      group = gp,
                      colour = src),
                  linewidth = 0.2) +
        scale_colour_manual(breaks = c("fb_d", "fb_r", "sim_d", "sim_r"),
                            labels = c("Data (seeder)",
                                       "Data (contact)",
                                       "Simulation (seeder)",
                                       "Simulation (contact)"),
                            values = c("red", "blue", "pink", "lightblue")) +
        lims(y = 0:1) +
        labs(x = "Time (days)",
             y = "Proportion",
             # title = "Kaplan-Meier plot by seeders / contact fish",
             title = str_glue("KM by seeder / contact{um}"),
             subtitle = str_glue("{params$dataset}/{params$label}: {description}"),
             colour = "Source") +
        theme_bw() +
        theme(panel.background = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              legend.position = "bottom") +
        facet_grid(cols = vars(variable),
                   rows = vars(trial),
                   scales = "free_x",
                   labeller = labeller(
                       variable = c(Tinf = "Proportion of family uninfected vs time",
                                    Tsym = "Proportion of family with no symptoms vs time",
                                    RP   = "Proportion of family surviving vs time"),
                       trial = c("1" = "Trial 1", "2" = "Trial 2")))
    plt
    
    # gfx_dir <- str_glue("{dataset}/gfx")
    #
    # ggsave(str_glue("{gfx_dir}/{dataset}-s{scen}-KM_donors-trials.pdf"),
    #        plot = plt, width = 12, height = 8)
    # ggsave(str_glue("{gfx_dir}/{dataset}-s{scen}-KM_grid-trials.png"),
    #        plot = plt, width = 12, height = 8, dpi = 600)
    # 
    # plt2 <- plt +
    #     facet_wrap(. ~ variable,
    #                scales = "free_x",
    #                labeller = labeller(
    #                    variable = c(Tsym = "Proportion of family uninfected vs time",
    #                                 Tsym = "Proportion of family with no symptoms vs time",
    #                                 RP   = "Proportion of family surviving vs time")))
    # 
    # ggsave(str_glue("{gfx_dir}/{dataset}-s{scen}-KM_grid.pdf"),
    #        plot = plt2, width = 12, height = 6)
    # ggsave(str_glue("{gfx_dir}/{dataset}-s{scen}-KM_grid.png"),
    #        plot = plt2, width = 12, height = 6, dpi = 600)
    # 
    # plt2
    
    list(plt = plt, data = data_t1)
}

