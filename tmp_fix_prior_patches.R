library(data.table)
library(stringr)
library(purrr)

files <- list.files("datasets/sim-test-inf/results", full.names = TRUE) |>
    str_sort(numeric = TRUE)

walk(files, \(f) {
    x <- readRDS(f)
    y <- with(x$params, readRDS(str_glue("datasets/{patch_dataset}/results/{patch_name}.rds")))
    
    priors <- x$params$priors
    setkey(priors, NULL)
    priors[, I := .I]
    
    priors2 <- merge(
        priors[, .SD, .SDcols = -"true_val"],
        y$params$priors[, .(parameter, true_val)],
        by = "parameter"
    ) |>
        setcolorder("use", after = "true_val") |>
        setorder(I)
    
    priors2[, I := NULL]
    
    x$params$priors <- priors2
    
    saveRDS(x, f)
})

