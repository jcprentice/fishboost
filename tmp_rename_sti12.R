{
    library(data.table)
    library(stringr)
    library(purrr)
}

ds <- "sim-test-inf2"

files <- list.files(str_glue("datasets/{ds}/results"),
                    full.names = TRUE) |>
    str_sort(numeric = TRUE)

walk(files, \(f) {
    x <- readRDS(f)
    x$params$dataset <- ds
    dirs <- names(x$params) |> str_subset("dir")

    walk(dirs, \(d) {
        x$params[[d]] <<- x$params[[d]] |>
            str_split_1("/") |>
            {\(x) {x[[2]] <- ds; x}}() |>
            str_flatten("/")
    })

    saveRDS(x, f)
})
