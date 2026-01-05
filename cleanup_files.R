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

