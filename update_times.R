library(stringr)
library(purrr)
library(lubridate)

ds <- "fb-qtest"

files <- list.files(str_glue("datasets/ds/data"), "bici") |>
    str_remove(".bici") |>
    str_sort(numeric = TRUE)

walk(files, \(f) {
    base <- str_glue("datasets/{ds}")
    f1 <- str_glue("{base}/data/{f}.bici")
    f2 <- str_glue("{base}/data/{f}-out/description.txt")
    f3 <- str_glue("{base}/results/{f}.rds")
    
    t1 <- readLines(f2)[[2]] |>
        str_remove("File written: ") |>
        ymd_hms()
    
    t2 <- ymd_hms(file.mtime(f1))
    
    t12 <- ceiling(as.numeric(t2 - t1) * 60 * 60)
    
    x <- readRDS(f3)
    x$time_start <- t1
    x$time_end <- t2
    x$time_taken <- list(tic = c(elapsed = 0),
                         toc = c(elapsed = t12))
    
    saveRDS(x, f3)
})
