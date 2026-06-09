{
    library(data.table)
    library(stringr)
    library(purrr)
}

files <- list.files("datasets/sim-test-inf2/results",
                    full.names = TRUE) |>
    str_sort(numeric = TRUE)
files

s1 <- files |> str_subset("scen-1-")
walk(s1, \(rf) {
    res <- readRDS(rf)
    res$parameter_estimates[str_detect(parameter, "_[GE]_"),
                            true_val := 0]
    res$params$priors[str_detect(parameter, "_[GE]_"),
                      true_val := 0]
    saveRDS(res, rf)
})


s2 <- files |> str_subset("scen-2-")
walk(s2, \(rf) {
    res <- readRDS(rf)
    res$parameter_estimates[str_detect(parameter, "_[GE]_.?i.?"),
                            true_val := 0]
    res$params$priors[str_detect(parameter, "_[GE]_.?i.?"),
                      true_val := 0]
    saveRDS(res, rf)
})

s4 <- files |> str_subset("scen-4-")
walk(s4, \(rf) {
    res <- readRDS(rf)
    res$parameter_estimates[str_detect(parameter, "r_[GE]_"),
                            true_val := 0]
    res$params$priors[str_detect(parameter, "r_[GE]_"),
                      true_val := 0]
    saveRDS(res, rf)
})
