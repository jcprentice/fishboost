cleanup_bici_files <- function(params) {
    if (file.exists(params$config)) {
        message("Removed old BICI script")
        file.remove(params$config)
    }
    
    if (dir.exists(params$output_dir)) {
        message("Removed old BICI files")
        cmd <- str_glue("rm -rf {params$output_dir}")
        system(cmd)
    }
}

