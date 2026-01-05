{
    library(data.table)
    library(purrr)
    library(stringr)
    
    source("rebuild_posteriors.R")
}

combine_replicates(dataset = "fb-final2") {
    data_dir <- str_glue("datasets/{dataset}/data")
    
    scen_reps <- data_dir |>
        list.files() |>
        str_remove_all("scen-|-out|-data.tsv|.xml") |>
        unique() |>
        str_sort(numeric = TRUE) |>
        tail(1) |>
        str_split_1("-") |>
        as.integer()
    
    scens <- seq_len(scen_reps[[1]])
    reps  <- seq_len(scen_reps[[2]])
    
    # Move over all trace files into scen-x-1
    walk(scens, \(scen) {
        # scen <- 1
        out_dir <- str_glue("{data_dir}/scen-{scen}-1-out")
        walk(reps, \(rep) {
            # rep <- 2
            from_dir <- str_glue("{data_dir}/scen-{scen}-{rep}-out")
            from <- list.files(from_dir, "trace_\\d") |>
                str_sort(numeric = TRUE)
            nfrom <- length(from)
            if (nfrom == 0) return()
            
            offset <- nfrom * (rep - 1L) - 1L
            to <- str_c("trace_", seq_len(nfrom) + offset, ".tsv")
            # print(str_glue("Rep {rep} : {from} -> {to}"))
            file.rename(str_glue("{from_dir}/{from}"),
                        str_glue("{out_dir}/{to}"))
        })
    })
    
    # Merge trace_combine.tsv and extended_trace_combine.tsv files
    walk(scens, \(scen) {
        # scen <- 1
        walk(c("", "extended_"), \(fp) {
            # fp <- ""
            tc_files <- list.files(data_dir,
                                   str_glue("^{fp}trace_combine"),
                                   recursive = TRUE,
                                   full.names = TRUE) |>
                str_subset(str_glue("scen-{scen}"))
            
            if (length(tc_files) <= 1) return()
            message(str_glue(" - Rebuilding {fp}trace_combine.tsv file for s{scen}"))
            
            foo <- map(tc_files, fread) |> rbindlist()
            foo[, state := seq(0L, .N - 1L)]
            fwrite(foo, tc_files[[1]], sep = "\t")
            file.remove(tc_files[-1])
        })
    })
    
    walk(scens, \(scen) {
        # scen <- 1
        out_dir <- str_glue("{data_dir}/scen-{scen}-1-out")
        # ebv_files <- str_glue("{data_dir}/scen-{scen}-{reps}-out/ebvs.csv")
        ebv_files <- list.files(data_dir, "ebvs", recursive = TRUE, full.names = TRUE) |>
            str_subset(str_glue("scen-{scen}"))
        if (length(ebv_files) <= 1) return()
        message(str_glue(" - Rebuilding 's{scen} / ebvs.csv'"))
        
        ebvs <- map(ebv_files, fread) |>
        rbindlist()
        ebvs2 <- ebvs[, map(.SD, mean), id]
        fwrite(ebvs2, str_glue("{out_dir}/ebvs.csv"))
        file.remove(ebv_files[-1])
    })
    
    walk(scens, ~ rebuild_sire_posteriors(dataset, str_glue("scen-{.x}-1-out")))
    
    # Delete remaining files
    walk(scens, \(scen) {
        walk(reps, \(rep) {
            if (rep == 1) return()
            message(str_glue(" - Deleting files for scen-{scen}-{rep} ..."))
            
            out_dir <- str_glue("{data_dir}/scen-{scen}-{rep}-out")
            files <- list.files(out_dir, full.names = TRUE, recursive = TRUE)
            if (length(files) > 0) file.remove(files)
            if (dir.exists(out_dir)) {
                message(str_glue("     - Deleting output files"))
                file.remove(out_dir)
            }
            files <- list.files(data_dir, str_glue("-{scen}-{rep}(-|\\.)"), full.names = TRUE)
            if (length(files) > 0) {
                message(str_glue("     - Deleting input files"))
                file.remove(files)
            }
        })
    })
}
