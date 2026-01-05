{
    library(data.table)
    library(stringr)
    library(purrr)
    library(ggplot2)
    library(ggcorrplot)
}

ebv_ranges <- function(dataset = "fb-simple", scen = 1) {
    # dataset <- "fb-simple"; scen <- 1
    
    message(str_glue("Getting ranges for '{dataset} / s{scen}' ..."))

    data_dir <- str_glue("datasets/{dataset}/data")
    out_dir <- str_glue("{data_dir}/scen-{scen}-1-out")
    res_dir <- str_glue("datasets/{dataset}/results")
    
    setup <- readRDS(str_glue("{res_dir}/scen-{scen}-1.rds"))$params$setup
    
    fstr <- switch(setup, "fb1" = "1", "fb2" = "2", "12")
    fb_data <-  readRDS(str_glue("fb_data/fb_data{fstr}.rds"))[, .(id, donor)]

    # Need to find what order SIRE put the BVs are in (should match the one
    # we'd get if we patched from the ebvs.csv file
    # EBVs.csv adds a few extra columns that we don't want

    f <- str_glue("{out_dir}/ebvs.csv")
    ebv_names <- if (file.exists(f)) {
        fread(f, nrows = 0) |> names() |> str_remove(" EBV") |>
            str_replace_all("s_" = "sus_", "i_" = "inf_", "t_" = "tol")
    } else {
        c("id", "inf_g", "inf_e", "sus_g", "sus_e", "tol_g",  "tol_e")
    }
    
    
    files <- list.files(str_glue("{out_dir}/states"),
                        full.names = TRUE)
    
    ebvs <- map(files, \(f) {
        line1 <- str_which(readLines(f), "Individual effects")
        ebvs <- fread(f, skip = line1)
        setnames(ebvs, ebv_names)
        ebvs[, donor := fb_data$donor]
    }) |> rbindlist(idcol = "file")
    
    ebvs[, `:=`(sus_p = sus_g + sus_e,
                inf_p = inf_g + inf_e,
                tol_p = tol_g + tol_e)]
    
    qs <- c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1)
    ebv_qs <- ebvs[, .(qsa = quantile(sus_g, qs),
                       qse = quantile(sus_e, qs),
                       qsp = quantile(sus_p, qs),
                       qia = quantile(inf_g, qs),
                       qie = quantile(inf_e, qs),
                       qip = quantile(inf_p, qs),
                       qta = quantile(tol_g, qs),
                       qte = quantile(tol_e, qs),
                       qtp = quantile(tol_p, qs)),
                   donor]
    
    plt <- ggcorrplot(cor(ebvs[, .(sus_g, sus_e, inf_g, inf_e, tol_g, tol_e)]), method = "circle")
    
    min_s <- ebvs[, .SD[which.min(sus_g)], file][, tail(sort(table(id)), 5)]
    max_s <- ebvs[, .SD[which.max(sus_g)], file][, tail(sort(table(id)), 5)]
    min_i <- ebvs[, .SD[which.min(inf_g)], file][, tail(sort(table(id)), 5)]
    max_i <- ebvs[, .SD[which.max(inf_g)], file][, tail(sort(table(id)), 5)]
    min_t <- ebvs[, .SD[which.min(tol_g)], file][, tail(sort(table(id)), 5)]
    max_t <- ebvs[, .SD[which.max(tol_g)], file][, tail(sort(table(id)), 5)]
    
    mget(c("min_s", "max_s", "min_i", "max_i", "min_t", "max_t", "ebv_qs", "plt"))
}

dataset <- "fb-simple"; scens <- 1:10

scens <- setNames(scens, str_c("s", scens))



ranges <- map(scens, ~ ebv_ranges(dataset, .x))

