{
    library(purrr)
    library(stringr)
    source("flatten_bici_states.R")
}

cmd_args <- commandArgs(trailingOnly = TRUE)
if (length(cmd_args) > 0) {
    dataset <- cmd_args[[1]]
    num <- as.integer(cmd_args[[2]])
} else {
    dataset <- "fb-test"
    num <- 1L
}

switch(dataset,
       "fb-test" = {scens <- 1:5; reps <- 1:5},
       "sim-test2" = {scens <- 1:5; reps <- 1:5},
       {scens <- 1:5; reps <- 1:5})

name <- expand.grid(reps, scens, "scen") |>
    rev() |> apply(1, str_flatten, "-")

flatten_bici_states(dataset, name[[num]], "inf")

