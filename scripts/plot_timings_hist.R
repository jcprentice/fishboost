{
    library(data.table)
    library(stringr)
    library(purrr)
    library(ggplot2)
    library(ggeasy)
    library(ggh4x)
}

timings_hist <- function(dataset = "sim-test1",
                         scens = 1:9,
                         use_means = FALSE,
                         apply_gammas = 0) {
    
    # dataset <- "fb-simple"; scens <- 1:10
    # dataset <- "fb-final-old"; scens <- 1L
    # dataset <- "sim-test1"; scens <- 1:9
    # use_means <- FALSE; apply_gammas <- 0
    
    trials <- 12
    
    um <- if (use_means) "_pm" else "_ps"
    ag <- if (apply_gammas > 0) str_glue("_gam{apply_gammas}") else ""
    
    f <- str_glue("datasets/{dataset}/meta/km_data{um}{ag}.rds")
    km_data <- readRDS(f)
    data <- map(scens, \(i) {
        x <- km_data[[i]]$data
        if ("Tinf_sire" %in% names(x)) {
            x[src == "fb", Tinf := Tinf_sire]
        }
        x1 <- melt(x, measure.vars = c("Tinf", "Tsym", "RP"),
                  variable.name = "var")
        x1[, src2 := str_c(trial, "_", donor)]
        x1
    })
    
    descriptions <- map_chr(scens, ~ km_data[[.x]]$params$description)
    
    breaks <- expand.grid(trial = 1:2, donor = 0:1) |>
        apply(1, str_flatten, "_")
    
    labels <- expand.grid(#src = c("Data", "Simulation"),
        trial = c("Trial 1", "Trial 2"),
        donor = c("Contact", "Seeder")) |>
        apply(1, str_flatten, " ") |>
        setNames(breaks)
    
    
    plts <- map(scens, \(i) {
        data_i <- data[[i]]
        description <- descriptions[[i]]
        params <- km_data[[i]]$params
        
        um <- if (use_means) "simulations from posterior means" else "sampled from posterior"
        ag <- if (apply_gammas > 0) str_glue(", gamma={apply_gammas}")
        
        subtitle <- str_flatten(c(with(params, str_glue("{dataset}/{label}: {description}")),
                                  "\n",
                                  str_flatten_comma(c(um, ag))))
        
        data1 <- data_i[(var == "Tinf" & value < 160) |
                            (var == "RP" & value < 50) |
                            (var == "Tsym" & value < 160)]
        
        ggplot(data1, aes(x = value)) +
            geom_histogram(aes(y = after_stat(density),
                               fill = src),
                           position = "identity",
                           alpha = 0.5,
                           binwidth = 1) +
            # geom_density(aes(y = after_stat(density),
            #                  fill = src),
            #              colour = NA,
            #              # position = "dodge",
            #              position = "identity",
            #              alpha = 0.5) +
            scale_fill_manual("Source",
                              breaks = c("fb", "sim"),
                              labels = c("Data", "Simulation"),
                              values = c("blue", "red")) +
            labs(x = "Time (days)",
                 y = "Density",
                 title = str_glue("Histogram of event times"),
                 subtitle = subtitle) +
            theme_bw() +
            theme(legend.position = "bottom") +
            easy_remove_gridlines() +
            facet_grid2(cols = vars(var),
                        rows = vars(src2),
                        scales = "free",
                        independent = "y",
                        labeller = labeller(
                            var = c(Tinf = "Time of infection",
                                    Tsym = "Time to first symptoms",
                                    RP   = "Period from symptoms to death"),
                            # src = c(fb = "Data", sim = "Simulation")))
                            src2 = labels))
    })
    
    
    gfx_dir <- str_glue("datasets/{dataset}/gfx/timings")
    if (!dir.exists(gfx_dir)) {
        dir.create(gfx_dir)
    }
    
    walk(scens, \(i) {
        plt_str <- str_glue("{gfx_dir}/{dataset}-s{i}-timings{um}{ag}.png")
        ggsave(plt_str, plts[[i]], width = 9, height = 6, unit = "in")
        message("plotting ", plt_str)
    })
    
    plts
}

{
    # dataset <- "fb-final"; scens <- 1:8
    # dataset <- "fb-donors"; scens <- 1:3
    dataset <- "fb-lp"; scens <- 1:12
    use_means <- FALSE
    apply_gammas <- 0
}

plts <- timings_hist(dataset = dataset,
                     scens = scens,
                     use_means = use_means,
                     apply_gammas = apply_gammas)
