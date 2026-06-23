{
    library(data.table)
    library(purrr)
    library(stringr)
}

files <- list.files("datasets", ".rds",
                    recursive = TRUE, full.names = TRUE) |>
    str_subset("r0_selection|meta", negate = TRUE) |>
    str_sort(numeric = TRUE)

modded_res_files <- c()

# Results files
res_files <- str_subset(files, "results")
walk(res_files, \(f) {
    x <- readRDS(f)
    touched <- FALSE

    if ("Tsym" %in% x$params$timings) {
        touched <- TRUE
        message(f, ": timings")
        x$params$timings <- x$params$timing |>
            str_replace_all("Tsym", "Tsign")
    }

    if ("Tsym" %in% x$params$pass_events) {
        touched <- TRUE
        message(f, ": pass_events")
        x$params$pass_events <- x$params$pass_events |>
            str_replace_all("Tsym", "Tsign")
    }

    if ("Tsym" %in% names(x$popn)) {
        touched <- TRUE
        message(f, ": popn")
        x$popn[, Tsign := Tsym]
        x$popn[, Tsym := NULL]
    }

    if (touched) {
        modded_res_files <<- c(modded_res_files, f)
        saveRDS(x, f)
    }

}, .progress = TRUE)


# Summary files
modded_sum_files <- c()
sum_files <- str_subset(files, "summary")
walk(sum_files, \(f) {
    x <- readRDS(f)
    touched <- FALSE

    if ("transitions" %in% names(x)) {
        touched <- TRUE
        x$transitions |>
            setnames("Tsym", "Tsign", skip_absent = TRUE)
    }

    if (touched) {
        modded_sum_files <<- c(modded_sum_files, f)
        saveRDS(x, file = f)
    }
}, .progress = TRUE)

