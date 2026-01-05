cleanup_sire_files <- function(params) {
    if (params$DEBUG == TRUE) {
        message("DEBUG is on! Leaving old files in place.")
        return()
    }

    if (str_detect(params$sire_version, "sire"))
        cleanup_sire_files(params)
    else
        cleanup_bici_files(params)
}


cleanup_sire_files <- function(params) {
    if (params$DEBUG) {
        message("DEBUG is on! Leaving old files in place.")
        return()
    }

    state_files <- list.files(params$states_dir, full.names = TRUE)
    if (length(state_files) > 0) {
        file.remove(state_files)
        message("Removed old SIRE state files")
    }
    
    out_files <- list.files(params$output_dir, full.names = TRUE) |>
        str_subset("txt|tsv|csv")
    if (length(out_files) > 0) {
        file.remove(out_files)
        message("Removed old SIRE output files")
    }
    
    in_files <- str_c(str_glue("{params$data_dir}/{params$name}"),
                        c(".xml", "-data.tsv", "-H.tsv"))
    in_flag <- FALSE
    walk(in_files, \(f) {
        if (file.exists(f)) {
            file.remove(f)
            in_flag <<- TRUE
        }
    })
    if (in_flag) message("Removed old SIRE input files")
}


cleanup_bici_files <- function(params) {
    config <- str_glue("{params$config}.bici")
    if (file.exists(config)) {
        message("Removed old BICI script")
        file.remove(config)
    }
    
    if (dir.exists(params$output_dir)) {
        message("Removed old BICI files")
        cmd <- str_glue("rm -rf {params$output_dir}")
        system(cmd)
    }
}

