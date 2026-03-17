{
    library(data.table)
    library(purrr)
    library(stringr)
    library(purrr)
    library(HDInterval)
    library(ggplot2)
}

plot_chain_ribbons <- function(dataset = "fb-test", scen = 1, use_hpdi = TRUE) {
    if (FALSE) {
        dataset <- "fb-test"; scen <- 1; use_hpdi <- TRUE
    }

    files <- list.files(str_glue("{dataset}/data/scen-{scen}-1-out"),
                        pattern = "trace",
                        full.names = TRUE) |>
        str_subset("combine", negate = TRUE) |>
        str_sort(numeric = TRUE)

    if (length(files) == 0) {
        message("No files found!")
        return(NULL)
    }

    x <- map(files, \(x) fread(x)) |>
        rbindlist(idcol = "chain")

    x[, str_subset(names(x), "state|Group|^G_|L_|Prior|Posterior|Number|log") := NULL]
    x[, names(.SD) := map(.SD, as.double)]

    xs <- x[, map(.SD, mean), chain] |>
        melt(id.var = "chain", value.name = "mean")

    if (use_hpdi) {
        xs_lo <- x[, map(.SD, \(x) hdi(x, credMass = 0.95)[["lower"]]), chain] |>
            melt(id.var = "chain", value.name = "min")
        xs_hi <- x[, map(.SD, \(x) hdi(x, credMass = 0.95)[["upper"]]), chain] |>
            melt(id.var = "chain", value.name = "max")
    } else {
        xs_lo <- x[, map(.SD, quantile, 0.025), chain] |>
            melt(id.var = "chain", value.name = "min")
        xs_hi <- x[, map(.SD, quantile, 0.975), chain] |>
            melt(id.var = "chain", value.name = "max")
    }
    xs[, `:=`(min = xs_lo$min,
              max = xs_hi$max)]


    plt <- ggplot(xs, aes(x = chain, colour = variable)) +
        geom_line(aes(y = mean)) +
        geom_point(aes(y = mean)) +
        geom_ribbon(aes(ymin = min, ymax = max, fill = variable),
                    alpha = 0.2) +
        expand_limits(y = 0) +
        labs(x = "Chain",
             y = "Value") +
        theme(legend.position = "none") +
        facet_wrap(. ~ variable,
                   scales = "free")

    plt
}

# dataset <- "fb-final"; scen <- 2
#
# plt <- plot_chain_ribbons(dataset, scen)
#
# plt_str <- str_glue("{dataset}/gfx/chains/{dataset}-s{scen}-chain_variables.pdf")
# ggsave(plt_str, plt, width = 16, height = 12)
#
