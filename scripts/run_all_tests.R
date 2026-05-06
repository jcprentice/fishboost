{
    library(stringr)
    library(data.table)
    library(purrr)

    source("fix_H_links.R")
    source("scripts/check_convergence.R")
    source("scripts/posterior.R")
    source("scripts/get_covariances.R")
    source("scripts/plot_ebvs.R")
    source("scripts/plot_ebvs2.R")
    source("scripts/pars_errorbars.R")
    source("scripts/pars_bias.R")
    source("plotting/plot_chains.R")
    source("plotting/plot_correlations.R")
    source("scripts/tornadoes.R")
    source("scripts/km_plots.R")
    source("scripts/model_fit_dev.R")
}

run_all_tests <- function(dataset = "fb-test",
                          km = FALSE) {
    if (FALSE) {
        dataset <- "fb-test"
        dataset <- "sim-test-inf2"
        km <- FALSE
    }

    message(str_glue("Generating all outputs for '{dataset}'"))

    fix_H_links(dataset)

    sr <- list.files(str_glue("datasets/{dataset}/results")) |>
        str_sort(numeric = TRUE) |>
        str_extract_all("\\d+") |>
        map(as.integer)
    scens <- map_int(sr, 1)
    reps  <- map_int(sr, 2)

    # Check convergence statistics
    out <- check_convergence(dataset)
    print(out$summary)

    # Overall errorbar plot of data
    pars_errorbars(dataset)

    if (str_detect(dataset, "fb")) {
        c(plot_chains, get_posterior, get_covariances, plot_ebvs, plot_correlations) |>
            walk(\(f) walk2(scens, reps, possibly(~ f(dataset, .x, .y))))
    }

    if (str_detect(dataset, "sim")) {
        walk2(scens, reps, possibly(~ plot_ebvs(dataset, .x, .y)))
        pars_bias(dataset)
        tornadoes(dataset)
    }

    if (km) {
        opts <- list(n_plots = 50, post = "sample")
        plotopts <- c("keep_small_groups", "extreme_sires", "drop_donors",
                      "mean", "ribbon", "fb_only", "t1", "t2")

        km_plots(dataset, scens, simulate_new_data = "no", opts = opts,
                 plotopts = plotopts[c(4, 5)])

        km_plots(dataset, scens, simulate_new_data = "no", opts = opts,
                 plotopts = plotopts[c(2, 4, 5)])

        km_plots(dataset, scens, simulate_new_data = "no", opts = opts,
                 plotopts = plotopts[0])

        # Don't want the extreme sires results for donor/trial KM plots
        list.files(str_glue("datasets/{dataset}/gfx/km_dt"), "_es",
                   full.names = TRUE) |>
            file.remove()

        fit <- model_fit_dev(dataset)
    } else {
        fit <- NULL
    }

    list(out = out, fit = fit)
}

if (FALSE) {
    run_all_tests("fb-test")
    run_all_tests("fb-qtest")
    run_all_tests("fb-qtest-1e5")
    run_all_tests("fb-qtest-5e5")
}
