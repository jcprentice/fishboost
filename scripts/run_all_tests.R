{
    library(stringr)
    library(data.table)
    library(purrr)

    source("scripts/check_convergence.R")
    source("scripts/pars_errorbars.R")
    source("scripts/posterior.R")
    source("scripts/plot_ebvs.R")
    source("plotting/plot_chains.R")
    source("plotting/plot_correlations.R")
    source("scripts/km_plots.R")
}

run_all_tests <- function(dataset = "fb-qtest",
                          include_km_plots = FALSE) {
    if (FALSE) {
        dataset <- "fb-test"
    }

    message(str_glue("Generating all outputs for '{dataset}'"))

    scen_reps <- list.files(str_glue("datasets/{dataset}/results")) |>
        str_sort(numeric = TRUE) |>
        str_extract_all("\\d+") |>
        map(as.integer)

    scens <- map_int(scen_reps, 1)
    reps  <- map_int(scen_reps, 2)

    # Check convergence statistics
    out <- check_convergence(dataset)
    print(out$summary)

    # Overall errorbar plot of data
    pars_errorbars(dataset)

    walk2(scens, reps, possibly(~ plot_chains(dataset, .x, .y)))
    walk2(scens, reps, possibly(~ get_posterior(dataset, .x, .y)))
    walk2(scens, reps, possibly(~ plot_ebvs(dataset, .x, .y)))
    walk2(scens, reps, possibly(~ plot_correlations(dataset, .x, .y)))

    if (include_km_plots) {
        f <- str_glue("datasets/{dataset}/meta/km_data_ps.rds")

        opts <- list(n_plots = 50,
                     use_means = FALSE)
        plotopts <- c("drop_small_groups",
                      "drop_donors",
                      "use_sire_Tinfs")[0]

        km_plots(dataset, scens,
                 simulate_new_data = "bici", opts, plotopts)
    }
    out
}

if (FALSE) {
    out <- c("fb-test", "fb-qtest") |>
        map(run_all_tests)
}
