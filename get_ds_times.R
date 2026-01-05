{
    library(data.table)
    library(stringr)
    library(purrr)
    library(lubridate)
}

my_seconds_to_time <- function(t) {
    as.character(seconds_to_period(t)) |>
        str_replace_all(c("H" = "h", "M" = "m", "S" = "s"))
}

get_ds_times <- function(dataset = "fb-final", mult = 1) {
    # dataset <- "fb-final"; mult <- 1
    
    files <- list.files(str_glue("datasets/{dataset}/results"), "scen",
                        full.names = TRUE) |>
        str_sort(numeric = TRUE)
    
    if (is_empty(files)) {
        message("No results files found!")
        return(NULL)
    }
    
    x <- map(files, \(f) {
        # "../scen-1-1.rds" -> "1-1"
        tmp <- str_split_i(f, "/", 4) |>
            str_remove_all("scen-|.rds") |>
            str_split_1("-") |>
            as.integer()
        Scen <- tmp[[1]]
        Rep <- tmp[[2]]
        
        if (!file.exists(f)) {
            return(list(Scenario = scen_str,
                        nsamples = NA_real_,
                        seconds = NA_integer_,
                        Time = NA_character_))
        }
        
        res <- readRDS(f)
        tt <- res$time_taken
        elapsed <- if ("toc" %in% names(tt)) {
            with(tt, unname(toc - tic))
        } else {
            tt[["elapsed"]]
        }
        t <- ceiling(mult * elapsed) |> as.integer()
        
        
        list(Scen = Scen, Rep = Rep,
             Samples = res$params$nsample,
             Seconds = t)
        # list(Scenario = scen_str,
        #      Samples = format(res$params$nsample,
        #                       scientific = FALSE,
        #                       big.mark = ","),
        #      Seconds = format(t,
        #                       scientific = FALSE,
        #                       big.mark = ","),
        #      Time = my_seconds_to_time(t))
    }) |>
        rbindlist()
    
    x[, .(Scenario = str_c(Scen, "-", Rep),
          Samples = format(Samples, scientific = FALSE, big.mark = ","),
          Seconds = format(Seconds, scientific = FALSE, big.mark = ","),
          Time = my_seconds_to_time(Seconds))] |>
        print()
    
    run_times <- map_dbl(files, ~ {
        tt <- readRDS(.x)$time_taken
        if ("toc" %in% names(tt)) tt$toc[["elapsed"]] else tt[["elapsed"]]
    })
    message(str_glue("Runs took about {a} (longest {b})",
                     a = seconds_to_period(ceiling(mean(run_times))),
                     b = seconds_to_period(ceiling(max(run_times)))))
    invisible(x)
}

# dataset <- "sim-test1"
# get_ds_times(dataset, mult = 1)

