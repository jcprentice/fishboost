{
    library(stringr)
    library(data.table)
    library(purrr)

    source("scripts/check_convergence.R")
    source("scripts/posterior.R")
    source("scripts/plot_ebvs.R")
    source("scripts/pars_errorbars.R")
    source("plotting/plot_chains.R")
    source("scripts/correlations.R")
    source("scripts/km_plots.R")
}

run_all_tests <- function(dataset = "fb-qtest",
                          include_km_plots = FALSE) {
    if (FALSE) {
        dataset <- "fb-qtest"
    }

    message(str_glue("Generating all outputs for '{dataset}'"))

    scens <- list.files(str_glue("datasets/{dataset}/results")) |>
        str_remove_all("scen-|.rds") |>
        str_sort(numeric = TRUE) |>
        str_split_i("-", 1) |>
        map_int(as.integer)

    reps <- list.files(str_glue("datasets/{dataset}/results")) |>
        str_remove_all("scen-|.rds") |>
        str_sort(numeric = TRUE) |>
        str_split_i("-", 2) |>
        map_int(as.integer)

    # Check convergence statistics
    out <- check_convergence(dataset)
    print(out$summary)

    # Overall errorbar plot of data
    pars_errorbars(dataset)

    walk2(scens, reps, {
        possibly(~ plot_chains(dataset, .x, .y))
        possibly(~ get_posterior(dataset, .x, .y))
        possibly(~ plot_ebvs(dataset, .x, .y))
        possibly(~ get_correlations(dataset, .x, .y))
    })

    if (include_km_plots) {
        f <- str_glue("datasets/{dataset}/meta/km_data_ps.rds")

        opts = list(n_plots = 50,
                    use_means = FALSE)
        plotopts = c("drop_small_groups",
                     "drop_donors",
                     "use_sire_Tinfs")[0]

        km_plots(dataset, scens,
                 simulate_new_data = "bici", opts, plotopts)
    }
    out
}

if (FALSE) {
    c("fb-test",
      "fb-qtest",
      "fb-qtest-1e5",
      "fb-qtest-5e5") |>
        walk(run_all_tests)
}
