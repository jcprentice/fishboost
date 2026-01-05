{
    library(data.table)
    library(purrr)
    library(stringr)
}

files <- list.files("datasets", recursive = TRUE, full.names = TRUE) |>
    str_subset("results.*\\.rds")

walk(files, \(f) {
    x <- readRDS(f)
    
    if ("sim_new_data" %notin% names(x$params)) {
        if ("use_fb_data" %notin% names(x$params)) {
            return()
        }
        
        message(str_glue("writing to '{f}'"))
        names(x$params)[[str_which(names(x$params), "use_fb_data")]] <- "sim_new_data"
        x$params$sim_new_data <- if (x$params$sim_new_data) "no" else "r"
        saveRDS(x, f)
        return()
        
    } else {
        y <- x$params$sim_new_data
        if (is.logical(y)) {
            message(str_glue("writing to '{f}'"))
            x$params$sim_new_data <- if (y) "r" else "no"
            saveRDS(x, f)
        }
    }
})
