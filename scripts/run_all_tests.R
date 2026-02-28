{
    library(stringr)
    library(data.table)
    library(purrr)

    source("fix_H_links.R")
    source("scripts/check_convergence.R")
    source("scripts/posterior.R")
    source("scripts/plot_ebvs.R")
    source("scripts/pars_errorbars.R")
    source("plotting/plot_chains.R")
    source("plotting/plot_correlations.R")
    source("scripts/tornadoes.R")
    source("scripts/km_plots.R")
    source("scripts/model_fit_rmsd.R")
}

run_all_tests <- function(dataset = "fb-qtest",
                          km = FALSE) {
    if (FALSE) {
        dataset <- "fb-test"
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
        walk2(scens, reps, possibly(~ plot_chains(dataset, .x, .y)))
        walk2(scens, reps, possibly(~ get_posterior(dataset, .x, .y)))
        walk2(scens, reps, possibly(~ plot_ebvs(dataset, .x, .y)))
        walk2(scens, reps, possibly(~ plot_correlations(dataset, .x, .y)))
    }

    if (str_detect(dataset, "sim")) {
        tornadoes(dataset)
    }

    if (km) {
        f <- str_glue("datasets/{dataset}/meta/km_data_ps.rds")

        opts = list(n_plots = 50,
                    use_means = FALSE)
        plotopts = c("drop_small_groups",
                     "drop_donors",
                     "use_sire_Tinfs")[0]

        km_plots(dataset, scens,
                 simulate_new_data = "bici", opts, plotopts)
        fit <- model_fit_rmsd(dataset)
    } else {
        fit <- NULL
    }

    list(out = out, fit = fit)
}

if (FALSE) {
    c("fb-test",
      "fb-qtest",
      "fb-qtest-1e5",
      "fb-qtest-5e5") |>
        walk(run_all_tests)
}
