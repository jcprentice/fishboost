source("plotting/tornadoes.R")

# run tests ----

# plts1 <- plot_tornadoes(data_set = "sim", scenario = 1, show_plots = FALSE)
# plts1 <- plot_tornadoes(data_set = "fb_1_12_Dx_linked", scenario = 17, show_plots = FALSE)

# data_set <- "fb-mpi"; combine <- FALSE
data_set <- "sim-censored-1"; combine <- FALSE
# data_set <- "sim-Gsi_cov_Da-1-mpi"; combine = FALSE
# data_set = "sim-donor_links1-2-mpi"; combine = TRUE

if (TRUE) {
    plts_sim <- list()
    for (i in 1:3) {
        # combine: 1=FALSE, 2=TRUE
        plts_sim[[i]] <- plot_tornadoes(data_set = data_set, combine = combine, scenario = i)
        
        # p_ht <- ceiling(sqrt(length(plts_sim[[i]]$plots))) * 5/3
        
        ggsave(glue("gfx/{data_set}/{data_set}-pars-scen{i}.pdf"),
               plts_sim[[i]]$pars,
               width = 10, height = 6)
    }
}

