source("plotting/posteriors.R")


# data_set <- "sim-censored-1"; num_scenarios <- 3
data_set <- "fb-parasites4"; num_scenarios <- 6

plts <- list()
for (i in 1:num_scenarios) {
    plts[[i]] <- plot_posteriors(data_set = data_set, scenario = i, ci = "hpdi")
    
    # n <- ceiling(length(plts[[i]]$plts) / 3)
    
    pdf_str = glue("gfx/{data_set}/{data_set}-posteriors-scen{i}.pdf")
    ggsave(pdf_str, plts[[i]]$pars, width = 10, height = 12)
    message("plotted ", pdf_str)
}

