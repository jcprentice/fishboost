library(data.table)
library(stringr)

export_trace <- function(dataset, scen) {
    # dataset <- "fb-final"; scen <- 2

    tc_txt <- str_glue("datasets/{dataset}/data/scen-{scen}-1-out/trace_combine.txt")
    tc <- str_replace(tc_txt, "\\.txt", "tsv")
    if (file.exists(tc_txt)) {
        message(" - renaming .txt to .tsv")
        file.rename(tc_txt, tc)
    }
    
    if (!file.exists(tc)) {
        stop(str_glue("{tc} missing!"))
    }
    
    x <- fread(tc)
    x[, str_subset(names(x), "Group") := NULL]
    fwrite(x, file = str_glue("{dataset}-{scen}.tsv"), sep = "\t")
}

export_trace("fb-final", 2)

