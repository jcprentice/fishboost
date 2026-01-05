{
    library(data.table)
    library(purrr)
    library(stringr)
    library(ggplot2)
    library(cowplot)
}

ind_traces <- function(dataset = "fb-fes3", scen = 1) {
    # dataset <- "fb-final-old"; scen <- 1
    message(str_glue("Getting values for {dataset} / s{scen}"))
    
    # name of state file
    sfiles <- list.files(str_glue("{dataset}/data/scen-{scen}-1-out"),
                         pattern = "state",
                         recursive = TRUE,
                         full.names = TRUE) |>
        str_sort(numeric = TRUE)
    
    event_times <- c("tE", "tI", "tD", "tR")
    
    map(seq_along(sfiles), \(i) {
        sfile <- sfiles[[i]]
        # search for line numbers featuring the word "Individual"
        lines <- str_which(readLines(sfile), "Individual")
        
        # fread only lines in between, fill columns
        x <- fread(sfile, skip = lines[[1]], nrows = diff(lines) - 1, fill = TRUE)
        setnames(x, c("id", event_times))
        x[, id := NULL]
        
        # Convert to numeric, with "no" as NA
        walk(event_times, \(e) {
            x[get(e) == "no", (e) := ""]
            x[, (e) := as.numeric(get(e))]
        })
        
        cbind(
            x[, map(.SD, max,  na.rm = TRUE)] |> setnames(str_c(event_times, "_max")),
            x[, map(.SD, mean, na.rm = TRUE)] |> setnames(str_c(event_times, "_mean"))
        )
    },
    .progress = TRUE) |>
        rbindlist()
}

dataset <- "fb-final"; nScens <- 8

tmp <- map(seq_len(nScens), \(i) ind_traces(dataset, i)) |>
    rbindlist(idcol = "scen")

maxes <- tmp[, .SD, .SDcols = str_subset(names(tmp), "scen|max")]
means <- tmp[, .SD, .SDcols = str_subset(names(tmp), "scen|mean")]
maxes[, id := seq_len(.N), scen]
means[, id := seq_len(.N), scen]

map_chr(seq_len(nScens), \(i) {
    f <- str_glue("datasets/{dataset}/results/scen-{i}-1.rds")
    title_str <- readRDS(f)$params$description |>
        str_remove(", convergence")
})



plt_tD <- map(seq_len(nScens), possibly(\(i) {
    f <- str_glue("datasets/{dataset}/results/scen-{i}-1.rds")
    params <- readRDS(f)$params
    title_str <- str_remove(params$description, ", convergence")
    ggplot() +
        geom_line(data = maxes, aes(x = id, y = tD_max), colour = "red") +
        geom_line(data = means, aes(x = id, y = tD_mean), colour = "blue") +
        expand_limits(y = 0) +
        labs(title = title_str)
}))

p <- plot_grid(plotlist = plt_tD, ncol = 2, byrow = FALSE)
ggsave(str_glue("datasets/{dataset}/gfx/tDs.png"),
       plot = p, width = 10, height = 10)

