{
    library(data.table)
    library(purrr)
    library(stringr)
}

get_subgroup_vals <- function(dataset = "fb-test", scen = 1) {
    if (FALSE) {
        dataset <- "fb-test"; scen <- 7
    }

    f <- str_glue("datasets/{dataset}/results/scen-{scen}-1.rds")
    popn <- readRDS(f)$popn
    pe <- readRDS(f)$parameter_estimates[str_starts(parameter, "beta|.P"), .(parameter, mean)]

    x <- popn[sdp == "progeny",
              .(id, trial, donor = fifelse(donor == 1, "Don", "Rec"))] |>
        _[, .N, .(trial, donor)]

    x[, `:=`(p = N / sum(N), N = NULL)] # should be 0.5, 0.5
    x[, `:=`(
        beta = pe[str_starts(parameter, "beta"), mean] |> rep(each = 2),
        LP   = pe[str_starts(parameter, "LP"), mean],
        DP   = pe[str_starts(parameter, "DP"), mean],
        RP   = pe[str_starts(parameter, "RP"), mean]
    )]

    # Summary of all subgroup values
    x[, .(trial, donor, beta = round(beta, 2),
          LP = round(LP, 2), DP = round(DP, 2), RP = round(RP, 2))] |>
        print()

    #
    x[, map(.SD, ~ sum(p * .x)) |> map(round, 2), .SDcols = -(1:3)] |>
        unlist()
}

get_subgroup_vals("fb-test", 1)
get_subgroup_vals("fb-test", 2)
get_subgroup_vals("fb-test", 7)
