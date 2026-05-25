{
    library(data.table)
    library(purrr)
    library(stringr)
}

files <- list.files("datasets", ".rds",
                    recursive = TRUE, full.names = TRUE) |>
    str_subset("r0_selection|meta", negate = TRUE) |>
    str_sort(numeric = TRUE)


# Results files
res_files <- str_subset(files, "results")
walk(res_files, \(f) {
    x <- readRDS(f)
    x$params$timings <- x$params$timing |>
        str_replace_all("Tsym", "Tsign")
    x$params$pass_events <- x$params$pass_events |>
        str_replace_all("Tsym", "Tsign")
    x$popn |>
        setnames("Tsym", "Tsign", skip_absent = TRUE)
}, .progress = TRUE)


# Summary files
sum_files <- str_subset(files, "summary")
walk(sum_files, \(f) {
    x <- readRDS(f)
    touched <- FALSE
    if ("transitions" %in% names(x)) {
        x$transitions |>
            setnames("Tsym", "Tsign", skip_absent = TRUE)
        touced <- TRUE
    }

    if (touched) saveRDS(x, file = f)
}, .progress = TRUE)

