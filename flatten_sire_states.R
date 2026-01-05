{
    library(data.table)
    library(purrr)
    library(stringr)
}

# Take the BICI state files and create an `etc.rds` file that contains all the
# parameters and traits in a data.table for easy access.

flatten_sire_states <- function(dataset = "fb-final", name = "scen-1-1", prog = FALSE) {
    # dataset <- "fb-final"; name <- "scen-1-1"; prog <- TRUE
    
    message(str_glue("Flattening SIRE state files for '{dataset}/{name}'"))
    
    base_dir   <- str_glue("datasets/{dataset}")
    out_dir    <- str_glue("base_dir/data/{name}-out")
    states_dir <- str_glue("{out_dir}/states")
    res_dir    <- str_glue("{base_dir}/results")
    
    files <-  list.files(states_dir, "state", full.names = TRUE) |>
        str_sort(numeric = TRUE)
    
    if (is_empty(files)) {
        message(" - No state files found!")
        return()
    } else {
        message(str_glue(" - found {length(files)} state files"))
    }
    
    rf <- str_glue("{res_dir}/{name}.rds")
    base_popn <- if (file.exists(rf)) {
        readRDS(rf)$popn[, .(id, sire, dam, sdp, trial, group, weight, donor)]
    }
    
    # Read first file to get positions of each state contained within
    lines <- readLines(files[[1]])
    
    # The top of the file contains the map of name:id. Transitions and IEs both
    # give `name`, where we want `id`, so just sub in `ids` where appropriate.
    headers <- c(str_which(lines, ":"), length(lines) + 1L)
    
    
    # This will read all the lines between headers n and n+1
    get_section <- function(lines, name) {
        n <- str_which(lines[headers], name)
        lines[seq(headers[[n]] + 1, headers[[n + 1]] - 1)]
    }
    
    parts <- c("parameters", "transitions", "ies")
    
    ie_names <- fread(str_glue("{out_dir}/ebvs.csv"), nrows = 0) |>
        names()
    
    out <- map(files, \(f) {
        # f <- files[[1]]
        lines1 <- readLines(f)
        
        parameters <- fread(text = get_section(lines1, "Parameters")) |>
            set_names("parameter", "value")
        
        transitions <- fread(text = get_section(lines1, "Individual transition"),
                             fill = TRUE) |>
            set_names("id", "Tinf", "Tinc", "Tsym", "Tdeath")
        
        transitions[, names(.SD) := map(.SD, \(x) {
            x <- as.character(x)
            fifelse(x %in% c("NA", "no"), NA_character_, x) |>
                as.numeric()
        }), .SDcols = -1]
        popn <- merge(base_popn, transitions, by = "id")
        
        ies <- fread(text = get_section(lines1, "Individual effects")) |>
            setnames(ie_names) |>
            setcolorder(c("id", "sus_g", "sus_e", "inf_g", "inf_e", "tol_g", "tol_e"),
                        skip_absent = TRUE)
        
        popn <- merge(popn, ies, by = "id")
        
        mget(parts)
    })
    
    etc <- map(parts, ~ map(out, .x) |> rbindlist(idcol = "state")) |>
        setNames(parts)
    
    f <- str_glue("{out_dir}/etc_inf.rds")
    message(str_glue(" - Saving to '{f}'"))
    saveRDS(etc, f)
    
    invisible(etc)
}

