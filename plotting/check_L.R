{
    library(data.table)
    library(stringr)
    library(purrr)
    library(ggplot2)
}

check_L <- function(dataset = "fb-final", scen = 1) {
    # dataset <- "fb-final"; scen <- 1
    
    files <- list.files(str_glue("datasets/{dataset}/data/scen-{scen}-1-out"),
                        pattern = "trace", full.names = TRUE) |>
        str_subset("combine", negate = TRUE) |>
        str_sort(numeric = TRUE)
    
    x <- map(files,
             \(f) fread(f)[, .(state,
                               L_inf_events,
                               L_trans_events)]) |>
        rbindlist(idcol = "chain")
    
    x[, state := state / state[2]]
    x1 <- x[state %% 100 == 0]
    x1[, `:=`(state = state / state[2],
              chain = as.factor(chain))]
    
    x2 <- melt(x1, measure.vars = c("L_inf_events", "L_trans_events"))
    
    plt <- ggplot(x2) +
        geom_line(aes(x = state, y = value, colour = chain)) +
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
 
