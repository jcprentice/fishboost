library(stringr)
library(purrr)
library(data.table)

dataset <- "fb-simple"

files <- list.files(dataset,
                    pattern = "RData",
                    recursive = TRUE,
                    full.names = TRUE)

walk(files, \(f) {
    RData_contents <- load(f)
    saveRDS(mget(RData_contents),
            file = str_replace(f, "RData", "rds"))
    file.remove(f)
})

