{
    library(stringr)
    library(purrr)
    source("run_bici_sim.R")
}

# Call Rscript "sim-base-inf" 1

cmd_args <- commandArgs(trailingOnly = TRUE)
run_from_script <- length(cmd_args) > 0

dataset <- if (run_from_script) cmd_args[[1L]] else "fb-test"

scens <- if (run_from_script) {
    as.integer(cmd_args[[2L]]) 
} else {
    list.files(str_glue("datasets/{dataset}/results")) |>
        str_split_i("-", 2) |>
        as.integer()
}

for (scen in scens) {
    reps <- list.files(str_glue("datasets/{dataset}/results"),
                       str_glue("scen-{scen}-")) |>
        str_remove_all("scen-.*-|.rds") |>
        as.integer() |>
        sort()
    
    for (rep in reps) {
        out <- run_bici_sim(dataset = dataset,
                            name = str_glue("scen-{scen}-{rep}"),
                            bici_cmd = "post-sim",
                            nreps = 50L)
        
        if (!is.null(out)) break
    }
}

# run_bici_sim("sim-base-inf", "scen-5-1", "post-sim", 50L)

