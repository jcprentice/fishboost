library(data.table)
library(stringr)

get_pars <- function(dataset = "testing", scen = 1) {
    f <- str_glue("datasets/{dataset}/results/scen-{scen}-1.rds")
    pe <- readRDS(f)$parameter_estimates
    
    pe[str_starts(parameter, "beta|latent|detection"),
       .(par = str_replace_all(parameter,
                               c("latent_period" = "LP",
                                 "detection_period" = "DP",
                                 "removal_period" = "RP")),
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
    x2[, `:=`(LP = pe[str_starts(parameter, "latent"), mean],
              DP = pe[str_starts(parameter, "detection"), mean],
              RP = pe[str_starts(parameter, "removal"), mean])]
    XPs <- x2[, .(LP = sum(LP * p) |> round(2),
                  DP = sum(DP * p) |> round(2),
                  RP = sum(RP * p) |> round(2))] |>
        as.list() |> unlist()
    
    c(mean_beta, XPs)
}

get_pars("testing", 1)
