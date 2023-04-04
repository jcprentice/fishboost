source("libraries.R")
source("source_files.R")

library(ggplot2)
library(cowplot)

# Generate data sets ----
generate_km_data <- function(data_set = "fb-parasites",
                             scenario = 3,
                             n_plots = 30,
                             apply_gammas = 5,
                             use_means = FALSE) {
    
    # Set to TRUE to use parameters while testing
    if (FALSE) {
        data_set = "fb-parasites4"
        scenario = 3
        n_plots = 3
        apply_gammas = 5
        use_means = FALSE
    }
    
    # Get parameters from target data set (and delete the bits we don't need)
    load(glue("results/{data_set}/scen-{scenario}-1.RData"))
    rm(list = c("pop", "parameter_estimates", "estimated_BVs", "pred_accs", "time_taken"))
    
    # Now override parameters
    params$use_fb_data     <- FALSE
    params$patch_data_set  <- data_set
    params$patch_scenario  <- scenario
    params$patch_with_mean <- use_means
    params$patch_traits    <- TRUE
    
    
    # Create empty list to store plots
    pop <- vector(mode = "list", length = n_plots + 1)
    
    for (i in 1:n_plots) {
        # This updates params with a different row from the trace file each time
        params <- patch_params(params)
        
        # Link traits (needs to be done after params is constructed)
        params <- apply_links(params)
        
        # Add in some Gamma distributed waiting times
        # Careful not to repeat this if not patching params
        if (apply_gammas > 0) {
            params$r_eta_rate  <- apply_gammas * params$r_eta_rate
            params$r_eta_shape <- apply_gammas * params$r_eta_shape
            params$r_rho_rate  <- apply_gammas * params$r_rho_rate
            params$r_rho_shape <- apply_gammas * params$r_rho_shape
            # params$r_gamma_rate  <- params$r_gamma_rate
            # params$r_gamma_shape <- params$r_gamma_shape
        }
        
        ## Generate pedigree and traits ----
        {
            pedigree <- make_pedigree(params)
            if (params$patch_traits) {
                traits <- patch_in_traits(pedigree, params)
            } else {
                traits <- make_traits_from_pedigree(pedigree, params)
            }
            
            # Set groups, trial, donors, and group effect
            traits <- set_groups(traits, params)
            
            # Set donor and trial effects
            traits <- apply_fixed_effects(traits, params)
        }
        
        # Simulate new data ----
        {
            message(glue("{i} / {n_plots}"))
            tmp <- simulate_epidemic(traits, params)
            pop[[i]] <- tmp[sdp == "progeny", .(id = i, sire, trial, group, Tsym, Trec, RP = Trec - Tsym, src = "sim")]
        }
    }
    
    message("Simulations complete, adding FB data")
    
    # Load the actual data, labelling as "fb"
    fb_data <- readRDS("fb_data/fb_data12.rds")
    
    pop[[n_plots + 1]] <- fb_data[sdp == "progeny",
                                  .(id = n_plots + 1, sire, trial, group, Tsym, Trec, RP = Trec - Tsym, src = "fb")]
    
    # Combine simulated and FB data
    data <- rbindlist(pop)
    
    # Combine parts
    list(data = data,
         params = params,
         opts = list(data_set = data_set,
                     scenario = scenario,
                     n_plots = n_plots,
                     apply_gammas = apply_gammas))
}



# Make Kaplan-Meier plots ----
plot_km_data <- function(data_list,
                         KM_measure = "Tsym",
                         KM_var = "sire",
                         trials = 12) {
    
    # Extract parts from data_list
    params       <- data_list$params
    data_set     <- data_list$opts$data_set
    scenario     <- data_list$opts$cenario
    apply_gammas <- data_list$opts$apply_gammas
    n_plots      <- data_list$opts$n_plots
    
    # Make sure to copy this otherwise we'll end up modifying it
    data         <- copy(data_list$data)
    
    # Rename the 2 parts we want to plot
    setnames(data, c(KM_var, KM_measure), c("var", "KM"))
    
    # Make a copy of the data and add in rows representing t=0
    data_t0 <- data[, .(KM = c(0, KM),
                        src = src[c(1, 1:.N)]),
                    .(id, var, trial)]
    
    # Rename src
    data_t0[, src := paste0(src, trial)]
    
    # Split by trial
    if (trials %in% 1:2) {
        data_t0 <- data_t0[trial == trials]
    }
    
    
    # Sort incorrectly places NA values at start, set them to Inf so they go to
    # the end, sort so they end at the back, then set them back to NA
    data_t0[is.na(KM), KM := Inf]
    setkey(data_t0, id, var, KM)
    data_t0[KM == Inf, KM := NA]
    
    # Survival curves by id and sire
    data_t0[, survival := 1 - seq(0, 1, length.out = .N), .(id, var)]
    
    # Cut off the last 10% since it tends to dominate... or just stop at 160 days
    # Tmax <- data_t0[!is.na(Trec), quantile(Trec, 0.9)[[1]]]
    Tmax <- data_t0[id == max(id), max(KM, na.rm = TRUE)]
    
    # Create a column to group by
    data_t0[, gp := factor(paste0(sprintf("%02d", id), "_", sprintf("%02d", var), "_", trial))]
    
    # String interpolation for the title
    wt_str <- if (apply_gammas) "gamma" else "exponential"
    prsts <- params$use_parasites
    if (prsts == "") prsts <- "not used"
    
    KMm_str <- switch(KM_measure,
                      "Tsym" = "Time to symptoms",
                      "Trec" = "Time to death",
                      "RP" = "Survival after symptoms",
                      "Measure")
    KMv_str <- switch(KM_var, "sire" = "family", "group" = "tank", "group")
    
    
    # Plot
    plt <- ggplot(data_t0) +
        geom_line(aes(x = KM, y = survival, group = gp, colour = src)) +
        coord_cartesian(xlim = c(0, 200)) +
        labs(x = "Time of death (days)",
             y = "Survival",
             title = glue("{KMm_str} by {KMv_str}"),
             # title = paste0("KM plot for FB vs simulated data, parasites ", prsts, ", (", wt_str, " waiting times)"),
             colour = "Source") +
        scale_color_manual(breaks = c("fb1", "fb2", "sim1", "sim2"),
                           labels = c("FB Trial 1", "FB Trial 2", "Sim Trial 1", "Sim Trial 2"),
                           values = c("blue", "red", "lightblue", "pink"))
    theme(panel.background = element_blank())
    
    
    
    list(plt = plt, data = data_t0)
}


# Play with functions ----
data_gam <- list()
out_Tsym <- list()
out_Trec <- list()
out_RP <- list()

for (i in 1:3) {
    data_set <- "fb-parasites4"
    n_plots <- 30
    scenario <- i #3 * i
    KM_var = "sire"
    trials <- 12
    
    data_gam[[i]] <- generate_km_data(data_set = data_set,
                                      scenario = scenario,
                                      n_plots = n_plots,
                                      apply_gammas = 5,
                                      use_means = FALSE)
    
    out_Tsym[[i]] <- plot_km_data(data_gam[[i]], KM_measure = "Tsym", KM_var = KM_var, trials = trials)
    out_Trec[[i]] <- plot_km_data(data_gam[[i]], KM_measure = "Trec", KM_var = KM_var, trials = trials)
    out_RP[[i]]   <- plot_km_data(data_gam[[i]], KM_measure = "RP",   KM_var = KM_var, trials = trials)
    
    plt <- plot_grid(out_Tsym[[i]]$plt, out_RP[[i]]$plt, nrow = 1)
    pdf_str <- glue("gfx/{data_set}/{data_set}-KM-scen{scenario}-{KM_var}-trial{trials}.pdf")
    ggsave(pdf_str, plt, width = 10, height = 4, unit = "in")
    png_str <- sub("pdf$", "png", pdf_str)
    ggsave(png_str, plt, width = 10, height = 4, unit = "in")
}

save(data_gam, out_Tsym, out_Trec, out_RP,
     file = glue("results/{data_set}/KM-scen{scenario}.RData"))

plot_grid(out_Tsym[[3]]$plt, out_RP[[3]]$plt)
