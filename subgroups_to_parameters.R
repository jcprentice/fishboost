library(data.table)
library(stringr)

get_pars <- function(dataset = "testing", scen = 1) {
    f <- str_glue("datasets/{dataset}/results/scen-{scen}-1.rds")
    pe <- readRDS(f)$parameter_estimates

    pe[str_starts(parameter, "beta|LP|DP"),
       .(par = parameter,
         mean = round(mean, 3))] |>
        with(setNames(mean, par)) |>
        print()

    x <- readRDS("fb_data/fb_12_drop71.rds")

    x1 <- x[sdp == "progeny", .N, trial]
    x1[, `:=`(p = N / sum(N), N = NULL)]
    x1[, beta := pe[str_starts(parameter, "beta"), mean]]
    mean_beta <- x1[, .(beta = sum(beta * p) |> round(2))] |>
        as.list() |> unlist()

    x2 <- x[sdp == "progeny", .N, .(donor, trial)]
    x2[, `:=`(p = N / sum(N), N = NULL)]
    x2[, `:=`(LP = rep(pe[str_starts(parameter, "LP"), mean], length = 4),
              DP = rep(pe[str_starts(parameter, "DP"), mean], length = 4),
              RP = rep(pe[str_starts(parameter, "RP"), mean], length = 4))]
    XPs <- x2[, .(LP = sum(LP * p) |> round(2),
                  DP = sum(DP * p) |> round(2),
                  RP = sum(RP * p) |> round(2))] |>
        as.list() |> unlist()

    c(mean_beta, XPs)
}

get_pars("fb-test", 2)

