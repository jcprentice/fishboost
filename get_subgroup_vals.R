{
    library(data.table)
    library(purrr)
    library(stringr)
}

get_subgroup_vals <- function(dataset = "fb-test", scen = 1) {
    # dataset <- "fb-test"; scen <- 1
    
    f <- str_glue("datasets/{dataset}/results/scen-{scen}-1.rds")
    popn <- readRDS(f)$popn
    pe <- readRDS(f)$parameter_estimates
    
    # pe[str_starts(parameter, "beta|latent|detection"),
    #    .(par = str_replace_all(parameter,
    #                            c("latent_period" = "LP",
    #                              "detection_period" = "DP",
    #                              "removal_period" = "RP")),
    #      mean = round(mean, 3))] |>
    #     with(setNames(mean, par)) |>
    #     print()
    # 
    # message("")
    
    x <- (popn[, .(id, sdp, trial, donor)]
          [sdp == "progeny", .N, .(trial, donor)]
          [, donor := fifelse(donor == 1, "Don", "Rec")])
    
    x[, `:=`(p = N / sum(N), N = NULL)] # should be 0.5, 0.5
    x[, `:=`(
        beta = rep(pe[str_starts(parameter, "beta"), mean], each = 2),
        LP = pe[str_starts(parameter, "latent"), mean],
        DP = pe[str_starts(parameter, "detection"), mean],
        RP = pe[str_starts(parameter, "removal"), mean]
    )]
    
    x[, .(trial, donor, beta = round(beta, 2),
          LP = round(LP, 2), DP = round(DP, 2), RP = round(RP, 2))] |>
        print()
    
    message("")
    
    means <- x[, map2(p, .SD, ~ sum(.x * .y)) |> map(round, 2) |> setNames(names(.SD)),
      .SDcols = c("beta", "LP", "DP", "RP")]
    
    means[1] |> as.list() |> unlist() |> setNames(names(means))
}

get_subgroup_vals("fb-test", 1)
