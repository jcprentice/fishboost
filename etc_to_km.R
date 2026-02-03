{
    library(stringr)
    library(purrr)
    library(data.table)
}

etc_to_km <- function(dataset = "fb-test", scens = 0, opts, bici_cmd = "post-sim") {
    base_dir <- str_glue("datasets/{dataset}")
    data_dir <- str_glue("{base_dir}/data")
    res_dir  <- str_glue("{base_dir}/results")
    meta_dir <- str_glue("{base_dir}/meta")
    
    files <- list.files(data_dir,
                        str_glue("etc_{bc}.rds",
                                 bc = switch(bici_cmd,
                                             sim = "sim", "ps")),
                        recursive = TRUE,
                        full.names = TRUE) |>
        str_sort(numeric = TRUE)
    
    if (is_empty(files)) {
        message("No etc files found")
        return(NULL)
    }
    
    if (all(scens != 0)) {
        files <- files |>
            keep(~ .x |> str_split_i("/", 4) |> str_split_i("-", 2) |>
                     as.integer() |> is.element(scens))
    }
    
    # Filter down to just 1 rep / scen
    dt <- data.table(files)
    dt[, `:=`(scen = files |> str_split_i("/", 4) |> str_split_i("-", 2) |> as.integer(),
              name  = files |> str_extract("(scen.*)-out", 1))]
    names <- dt[, name[[1]], scen][, V1]
    
    km_data <- map(names, \(name) {
        # name <- names[[1]]
        ef <- str_glue("{data_dir}/{name}-out/etc_ps.rds")
        rf <- str_glue("{res_dir}/{name}.rds")
        
        if (!file.exists(ef) || !file.exists(rf)) return(NULL)
        
        etc <- readRDS(ef)$popn[sdp == "progeny",
                                .(id = state, sire, trial, donor, group,
                                  Tinf, Tsym, Tdeath, RP = Tdeath - Tsym,
                                  src = "sim")]
        n <- last(etc$id)
        res <- readRDS(rf)
        
        data <- rbind(etc,
                      res$popn[sdp == "progeny",
                               .(id = n + 1L, sire, trial, donor, group,
                                 Tinf, Tsym, Tdeath, RP = Tdeath - Tsym,
                                 src = "fb")])
        
        list(data = data, params = res$params, opts = opts)
    })
    
    if (!dir.exists(meta_dir)) {
        message("- mkdir ", meta_dir)
        dir.create(meta_dir)
    }
    
    saveRDS(km_data, file = str_glue("{meta_dir}/km_data_ps.rds"))
    
    km_data
}
