{
    library(purrr)
    library(stringr)

    source("flatten_bici_states.R")
}

cmd_args <- commandArgs(trailingOnly = TRUE)
n <- as.integer(cmd_args[[1]])

list.files("datasets/sim-test-inf1/results",
           full.names = TRUE) |>
    str_sort(numeric = TRUE) |>
    pluck(n) |>
    map(readRDS) |>
    map("params") |>
    walk(flatten_bici_states)
