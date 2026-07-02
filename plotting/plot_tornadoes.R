{
    library(data.table)
    library(purrr)
    library(stringr)
    library(ggplot2)
    library(ggtext)
    library(ggeasy)
    library(cowplot)

    source("rename_pars.R")
    source("utils.R")
}

plot_tornadoes <- function(dataset = "sim-test-inf1",
                           scen = 2) {

    if (FALSE) {
        dataset <- "sim-test-inf1"
        scen <- 2
    }

    message(str_glue("Generating tornado plots for {dataset} / s{scen} ..."))

    # Directories and list of files to scan
    data_dir    <- str_glue("datasets/{dataset}/data")
    results_dir <- str_glue("datasets/{dataset}/results")
    res_files   <- list.files(results_dir,
                              pattern = str_glue("scen-{scen}-"),
                              full.names = TRUE) |>
        str_sort(numeric = TRUE)

    x <- res_files |>
        map(readRDS) |>
        map("parameter_estimates") |>
        rbindlist(idcol = "id") |>
        _[!str_starts(parameter, "G_"),
           .(id, parameter, true_val, median, hdi95min, hdi95max)]

    pars <- unique(x$parameter)
    html_pars <- html_names(pars) |> setNames(pars)

    params <- readRDS(res_files[[1]])$params
    priors <- params$priors


    plts <- map(pars, \(par) {
        if (FALSE) {
            i <- 1; par <- pars[[i]]
        }

        xp <- x[parameter == par, .(mu = median, true_val, hdi95min, hdi95max)] |>
            setorder(mu)
        xp[, id := .I]

        xlines <- c(xp$true_val[[1]], mean(xp$mu),
                    hdi(xp$mu)[c("lower", "upper")])

        prior <- priors[parameter == par, .(val1, val2)] |> as.numeric()

        ggplot(xp) +
            geom_segment(aes(x = hdi95min, xend = hdi95max,
                             y = id),
                         colour = "red") +
            geom_point(aes(x = mu,
                           y = id),
                       colour = "red") +
            geom_vline(xintercept = xlines,
                       colour = c("green", "blue", "blue", "blue"),
                       linewidth = c(1, 1, 0.5, 0.5),
                       linetype = c("solid", "solid", "dashed", "dashed"),
                       show.legend = FALSE) +
            scale_x_continuous(limits = ~ range(.x, prior, 0)) +
            labs(x = "Value",
                 title = html_pars[[par]]) +
            theme_classic() +
            theme(plot.title = element_markdown()) +
            easy_remove_y_axis()
    }) |>
        setNames(pars)

    plts$title_plt <- ggplot() +
        labs(title = str_glue("{dataset} / s{scen}: {params$description}")) +
        theme_classic()

    plts$empty <- ggplot() + theme_classic()

    plts
}

