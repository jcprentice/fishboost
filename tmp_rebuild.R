{
    library(purrr)
    library(stringr)
    source("rebuild_posteriors.R")
}

datasets <- c("fb-test",
              "sim-events",
              "sim-test-inf1",
              "sim-test-inf2")

walk(datasets, \(ds) {
    names <- list.files(str_glue("datasets/{ds}/results")) |>
        str_remove_all(".rds") |>
        str_sort(numeric = TRUE)

    walk(names, ~ rebuild_bici_posteriors(ds, .x))
})
