library(data.table)
library(purrr)
library(stringr)

source("read_sire_state_file.R")

states_to_etc <- function(dataset = "fb-final", scens = 0) {
    # dataset <- "fb-final"; scens <- 0
    data_dir <- str_glue("datasets/{dataset}/data")
    
    if (length(scens) == 1 && scens == 0) {
        scens <- list.files(data_dir) |>
            str_remove_all("scen-|-out|.xml|-data.tsv") |>
            unique() |>
            str_split_i("-", 1) |>
            as.integer() |>
            sort()
    }
    
    walk(scens, \(scen) {
        # scen <- 1
        message("Rebuilding etc file for s", scen)
        out_dir <- str_glue("{data_dir}/scen-{scen}-1-out")
        
        etc_fn <- str_glue("{out_dir}/extended_trace_combine.tsv")
        if (file.exists(etc_fn)) return()
        
        
        states <- list.files(str_glue("{out_dir}/states"), full.names = TRUE)
        
        out <- map(states, \(state) {
            # state <- states[[1]]
            tmp <- read_sire_state_file(state, "parameters,ies")
            ies <- tmp$ies |> melt(id.vars = "id")
            ies[, `:=`(variable = str_c(id, "_", variable), id = NULL)]
            rbind(tmp$parameters, ies)
        }, .progress = TRUE) |>
            rbindlist(idcol = "state") |>
            dcast(state ~ variable)
        
        out[, state := seq(0, .N - 1)]
        
        fwrite(out, file = etc_fn, sep = "\t")
    })
}

states_to_etc(dataset = "fb-final", scens = 1:8)
