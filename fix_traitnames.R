{
    library(stringr)
    library(purrr)
    library(data.table)
}

files <- list.files("datasets", ".rds", recursive = TRUE, full.names = TRUE) |>
    str_subset("meta|r0|km|KM|etc", negate = TRUE) |>
    str_sort(numeric = TRUE)

traits <- c("sus", "inf", "lat", "det", "tol")
mats <- c("Sigma_G", "Sigma_E", "Sigma_P",
          "cor_G", "cor_E", "cor_P",
          "cov_G", "cov_E", "cov_P",
          "h2", "fe_vals")

rep_traits <- function(x) {
    x |> str_replace_all(c("susceptibility" = "sus",
                           "infectivity"    = "inf",
                           "latency"        = "lat",
                           "detectability"  = "det",
                           "tolerance"      = "tol"))
}

walk(files, \(f) {
    x <- readRDS(f)
    
    setnames(x$popn, rep_traits(names(x$popn)))
    
    walk(mats, \(mat) {
        if (mat %in% names(x$params)) {
            dimnames(x$params[[mat]]) <<- map(dimnames(x$params[[mat]]), rep_traits)
        }
    })
    
    saveRDS(x, f)
})

