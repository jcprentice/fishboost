{
    library(data.table)
    library(stringr)
    library(purrr)
    library(cowplot)

    source("utils.R")
    source("get_plot_matrix.R")
    source("plotting/plot_posteriors.R")
}

get_posterior <- function(dataset = "fb-test", scen = 1, rep = 1,
                          output = "png") {
    if (FALSE) {
        dataset <- "fb-test"; scen <- 1; rep <- 1
        dataset <- "sim-base-inf"; scen <- 1; rep <- 1
        output <- "png"
    }

    {
        base_dir <- str_glue("datasets/{dataset}")
        res_dir  <- str_glue("{base_dir}/results")
        gfx_dir  <- str_glue("{base_dir}/gfx")
        post_dir <- str_glue("{gfx_dir}/posteriors")
        cov_dir  <- str_glue("{gfx_dir}/posteriors_covs")

        c(gfx_dir, post_dir, cov_dir) |>
            discard(dir.exists) |>
            walk(~ message("- mkdir ", .x)) |>
            walk(dir.create, recursive = TRUE)

    }

    plts <- plot_posteriors(dataset, scen, rep, ci = "hpdi", draw = "density")

    pars <- names(plts$plts)
    pmat <- get_plot_matrix(pars, "traits")

    plt <- plot_grid(plts$title_plt,
                     plot_grid(plotlist = plts$plts[pmat$plt_names],
                               nrow = pmat$nr,
                               ncol = pmat$nc),
                     ncol = 1, rel_heights = c(0.3, pmat$nr))

    plt_str <- str_glue("{post_dir}/{dataset}-s{scen}-{rep}-posteriors")
    walk(output, ~ ggsave(str_c(plt_str, ".", .x), plt,
                          width = 2 * pmat$nc,
                          height = 2 * (pmat$nr + 0.3)))
    message(str_glue("- plotted '{plt_str}'"))


    # Covariance only
    cov_pars <- pars |> str_subset("_[GEP]_|h2_")

    if (length(cov_pars) == 0) {
        message("- No Genetic Variance parameters detected")
        return()
    }

    cpmat <- get_plot_matrix(cov_pars, "rect")

    cov_plt <- plot_grid(plts$title_plt,
                         plot_grid(plotlist = plts$plts[cpmat$plt_names],
                                   nrow = cpmat$nr,
                                   ncol = cpmat$nc),
                         ncol = 1, rel_heights = c(0.3, cpmat$nr))

    plt_str <- str_glue("{cov_dir}/{dataset}-s{scen}-{rep}-covs") |>
        str_replace("scen-", "s")
    walk(output, ~ ggsave(str_c(plt_str, ".", .x), cov_plt,
                          width = 2 * cpmat$nc,
                          height = 2 * (cpmat$nr + 0.3)))
    message(str_glue("- Plotted '{plt_str}'"))

    plt
}

if (FALSE) {
    # get_posteriors(dataset = "testing", scens = 1:2)
    dataset <- "fb-test-1e7"; scens <- 1:6
    walk(scens, possibly(\(i) get_posterior(dataset, i, 1, "png")))

    dataset <- "fb-qtest"; scens <- 1:18
    walk(scens, possibly(\(i) get_posterior(dataset, i, 1)))

    # get_posteriors(dataset = "sim-base-inf", scens = 5, reps = 1:5, combine = FALSE)
    # walk(scens, \(i) get_posterior(dataset, i, 1))
}
