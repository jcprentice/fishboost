library(stringr)
library(purrr)

datasets <- list.files("datasets")

walk(datasets, \(ds) {
    files <- list.files(str_glue("datasets/{ds}/results"), full.names = TRUE)

    walk(files, \(f) {
        x <- readRDS(f)
        names(x$params)[which(names(x$params) == "use_traits")] <- "ies"
        # saveRDS(x, f)
    })
})

