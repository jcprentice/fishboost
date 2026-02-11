{
    library(data.table)
    library(stringr)
    library(purrr)
    library(ggplot2)
}

check_L <- function(dataset = "fb-final", scen = 1) {
    if (FALSE) {
        dataset <- "fb-test"; scen <- 3
    }

    files <- list.files(str_glue("datasets/{dataset}/data/scen-{scen}-1-out/output-inf"),
                        pattern = "param_", full.names = TRUE) |>
        str_subset("combine", negate = TRUE) |>
        str_sort(numeric = TRUE)

    x <- map(files, fread) |>
        rbindlist(idcol = "chain")
    x[, str_subset(names(x), "chain|State|^L\\^", negate = TRUE) := NULL]
    x[, State := State / State[2]]

    # Thin down to at most 1e4 samples
    x1 <- x[seq(1, .N, length.out = 1e4) |> ceiling() |> unique()]
    x1[, `:=`(State = seq(.N),
              chain = as.factor(chain))]

    Lcols <- c("L^markov", "L^ie", "L^dist")
    x1[, names(.SD) := map(.SD, as.numeric), .SDcols = Lcols]
    x2 <- melt(x1, measure.vars = Lcols)

    plt <- ggplot(x2) +
        geom_line(aes(x = State, y = value, colour = chain)) +
        facet_wrap(. ~ variable, scales = "free_y", ncol = 1)

    plt
}

dataset = "fb-fes1"; scens = 1:10
# dataset = "fb-fes2"; scens = 1:11

plts <- map(n_scens, \(i) {
    plt <- check_L(dataset, i)
    ggsave(str_glue("datasets/{dataset}/gfx/L-scen{i}.png"),
           plot = plt, width = 12, height = 6)
    plt
})

