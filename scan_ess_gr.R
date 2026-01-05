{
    library(data.table)
    library(stringr)
    library(purrr)
}

scan_ess_gr <- function(dataset = "fb-test") {
    
    files <- list.files(str_glue("datasets/{dataset}/results"), "scen") |>
        str_remove("scen-") |>
        str_remove("\\.rds") |>
        str_sort(numeric = TRUE)
    
    out <- map(files, \(f) {
        f2 <- str_glue("datasets/{dataset}/results/scen-{f}.rds")
        if (!file.exists(f2)) return(data.table(ESS = NA, GR = NA))
        x <- readRDS(f2)
        x$parameter_estimates[!str_starts(parameter, "Group"),
                              .(scen = str_remove(x$params$name, "cen-"),
                                ESS = min(ESS, na.rm = TRUE),
                                GR = round(max(GR, na.rm = TRUE), 2))]
    }) |> rbindlist()
    print(out)
    
    out[, scen := files]
    out[, .(min_ESS = min(ESS), max_GR = max(GR))]
}

scan_ess_gr("fb-test")
scan_ess_gr("fb-dp")
scan_ess_gr("fb-donors")
scan_ess_gr("sim-test")
