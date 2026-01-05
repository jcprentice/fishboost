{
    library(stringr)
    library(purrr)
    source("run_bici_sim.R")
}

# Generate data sets ----
generate_km_data_bici <- function(dataset = "fb-test",
                                  scen = 1,
                                  opts = list(n_plots = 50L,
                                              use_means = FALSE)) {
    if (FALSE) {
        dataset <- "sim-base-inf"
        scen <- 5
        opts <- list(n_plots = 50L, use_means = FALSE)
    }
    
    data_dir <- str_glue("datasets/{dataset}/data")
    res_dir <- str_glue("datasets/{dataset}/results")
    
    n_plots <- opts$n_plots
    use_means <- opts$use_means
    
    f <- list.files(data_dir, "etc_ps.rds",
                    recursive = TRUE, full.names = TRUE) |>
        keep(~ .x |> str_split_i("/", 4) |> str_split_i("-", 2) |>
                 as.integer() |> is.element(scen))
    
    bici_data <- if (is_empty(f)) {
        files <- list.files(data_dir, str_glue("scen-{scen}-.*bici")) |>
            str_remove_all("\\.bici") |>
            str_sort(numeric = TRUE)
        
        for (f in files) {
            out <- run_bici_sim(dataset = dataset, name = f,
                                bici_cmd = "post-sim",
                                nreps = 0)
            if (!is.null(out)) break
        }
        out$popn
    } else {
        name <- str_split_i(f[[1]], "/", 4) |>
            str_remove("-out")
        readRDS(f[[1]])$popn
    }
    
    bici_data[, `:=`(id = state, src = "sim", state = NULL)]
    
    rf <- readRDS(str_glue("{res_dir}/{name}.rds"))
    rf$popn[, `:=`(id = n_plots + 1L,
                   src = "fb")]
    
    data <- rbind(bici_data,
                  rf$popn[, .SD, .SDcols = intersect(names(bici_data), names(rf$popn))],
                  fill = TRUE)
    
    data <- data[sdp == "progeny"]
    data[, sdp := NULL]
    data[, RP := Tdeath - Tsym]
    setcolorder(data, "RP", after = "Tdeath")
    
    list(data = data,
         params = rf$params,
         opts = c(dataset = dataset,
                  scen = scen, rep = rep,
                  opts))
}
