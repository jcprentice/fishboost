{
    library(data.table)
    library(stringr)
    library(purrr)
}

# Check convergence ----

check_convergence <- function(dataset = "fb-test") {
    # dataset <- "sim-base-inf"
    # dataset <- "fb-test"
    
    files <- list.files(str_glue("datasets/{dataset}/results"),
                        full.names = TRUE) |>
        str_sort(numeric = TRUE)
    
    pes <- map(files, readRDS) |>
        map("parameter_estimates") |>
        rbindlist(idcol = "file")
    
    {
        # Remove group effects
        pes <- pes[!str_starts(parameter, "G_")]
        
        # Add sce, rep, and description
        scens <- files |> str_split_i("/", 4) |> str_split_i("-", 2) |> as.integer()
        reps <- files |> str_split_i("/", 4) |> str_split_i("-", 3) |>
            str_remove_all(".rds") |> as.integer()
        descriptions <- map(files, readRDS) |> map("params") |> map_chr("description") |>
            str_remove_all("FB_12_rpw, |, GRM \\w*|, convergence")
        
        pes[, `:=`(scen = scens[file],
                   rep = reps[file],
                   desc = descriptions[file])] |>
            setcolorder(c("scen", "rep", "desc"))
        
        # Rename LPs
        pes[, parameter := parameter |>
                str_replace_all(c("latent_period" = "LP",
                                  "detection_period" = "DP",
                                  "removal_period" = "RP"))]
        
        pes[, GEV := desc |> str_split_1(", ") |> str_subset("GEV ") |> str_remove("GEV "), .I]
        pes[is.na(GEV), GEV := "SIT"]
        pes[, FE := desc |> str_split_1(", ") |> str_subset("FE  ") |> str_remove("FE "), .I]
        pes[is.na(FE), FE := "SIT"]
        
        # Tidy up
        pes[, `:=`(file = NULL, ci95min = NULL, ci95max = NULL)]
        setcolorder(pes, c("scen", "rep", "GEV"))
        }
    
    pes2 <- pes[, .(desc = first(desc),
                    # GEV = first(GEV),
                    # FE = first(FE),
                    mean_ESS = mean(ESS, na.rm = TRUE) |> round(),
                    min_ESS = ESS |> min(na.rm = TRUE),
                    mean_GR = mean(GR, na.rm = TRUE) |> round(2),
                    max_GR = GR |> max(na.rm = TRUE) |> round(2)),
                .(scen, rep)]
    pes2[, converged := fcase(min_ESS >= 500 & max_GR < 1.05, "***",
                              min_ESS >= 200 & max_GR < 1.1,  "**",
                              min_ESS >= 100 & max_GR < 1.2,  "*",
                              default = "")]
    
    list(summary = pes2,
         worst = pes[order(ESS)[1:20]])
}

