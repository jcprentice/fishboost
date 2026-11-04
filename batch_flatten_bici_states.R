{
    library(purrr)
    library(stringr)

    safe_source("flatten_bici_states.R")
}

cmd_args <- commandArgs(trailingOnly = TRUE)
if (length(cmd_args) > 0) {
    dataset <- cmd_args[[1]]
    num <- as.integer(cmd_args[[2]])
} else {
    dataset <- "sim-events"
    num <- 1L
}

name <- list.files(str_glue("datasets/{dataset}/data")) |>
    str_subset("bici") |>
    str_remove("\\.bici") |>
    str_sort(numeric = TRUE) |>
    _[[num]]

message(str_glue("Flattening {dataset}, {name}"))

params <- list(dataset = dataset,
               name = name,
               bici_cmd = "inf")

flatten_bici_states(params)

message("Finished!")
