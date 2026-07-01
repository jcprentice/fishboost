{
    library(purrr)
    library(stringr)

    source("plotting/plot_tornadoes.R")
    source("utils.R")
    source("get_plot_matrix.R")
}

tornadoes <- function(dataset = "fb-test", scens = 0, combine = TRUE) {
    if (FALSE) {
        dataset <- "sim-test-inf1"
        scens <- 0
        combine <- TRUE
    }

    gfx_dir <- str_glue("datasets/{dataset}/gfx/tornadoes")
    if (!dir.exists(gfx_dir)) {
        message("- mkdir ", gfx_dir)
        dir.create(gfx_dir, recursive = TRUE)
    }

    if (all(scens == 0)) {
        scens <- list.files(str_glue("datasets/{dataset}/results")) |>
            str_sort(numeric = TRUE) |>
            str_extract("-(\\d+)-", group = 1) |>
            unique() |>
            as.integer()
    }

    plts <- map(scens, \(scen)
                plot_tornadoes(dataset = dataset,
                               scen = scen)) |>
        setNames(str_c("s", scens))


    # Create a plot with subplots aligned by trait

    walk2(plts, names(plts), \(x, scen) {
        if (FALSE) {
            scen <- "s1"
            x <- plts[[scen]]
        }
        if (is.null(x)) return()

        pmat <- get_plot_matrix(names(x), "traits")

        plt <- plot_grid(x$title_plt,
                         plot_grid(plotlist = x[pmat$plt_names],
                                   nrow = pmat$nr,
                                   ncol = pmat$nc),
                         ncol = 1,
                         rel_heights = c(0.3, pmat$nr))

        plt_str <- str_glue("{gfx_dir}/{dataset}-{scen}-tornadoes.pdf")
        ggsave(plt_str, plt,
               width = 2.4 * pmat$nc,
               height = 2 * (pmat$nr + 0.3))
        message(str_glue("plotted '{plt_str}'"))
    })
}

# run tests ----
if (FALSE) {
    # tornandoes(dataset = "sim-test2", scens = 1:5, combine = TRUE)
    tornadoes("sim-base-inf")
    tornadoes("sim-test-inf1")
    tornadoes("sim-test-inf2")
}
