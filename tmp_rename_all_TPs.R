{
    library(data.table)
    library(stringr)
    library(purrr)
}

dataset <- "fb-qtest"

files <- list.files(str_glue("datasets/{dataset}/results")) |>
    str_sort(numeric = TRUE)

rename_TPs <- function(x) {
    x |> str_replace_all(c(
        "latent_period_"    = "LP_",
        "detection_period_" = "DP_",
        "removal_period_"   = "RP_"
    ))
}

walk(files, \(f) {
    x <- readRDS(f)

    x$params$priors[, parameter := rename_TPs(parameter)]
    x$parameter_estimates[, parameter := rename_TPs(parameter)]

    saveRDS(x, f)
})
