{
    library(stringr)
    library(purrr)
}

fix_H_links <- function(...) {
    pln <- if (Sys.info()[["sysname"]] == "Darwin") "gln" else "ln"
    
    walk(list(...), \(d) {
        base_dir <- str_glue("datasets/{d}/data")
        files <- list.files(base_dir, "^H", recursive = TRUE, full.names = TRUE)
        
        walk(files, \(f) {
            src <- str_glue("fb_data/{mat}", mat = str_split_i(f, "/", 5))
            if (!file.exists(src)) {
                message(str_glue("missing {src}"))
                return()
            }
            cmd <- str_glue("{pln} -sfr {src} {f}")
            system(cmd)
        })
    })
}

