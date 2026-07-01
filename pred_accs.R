{
    library(data.table)
    library(stringr)
    library(purrr)
    library(ggplot2)
    library(cowplot)
    library(pipebind)
}

pred_accs_plot <- function(dataset = "sim-test-inf", scen = 1, rep = 1, parents_only = TRUE) {
    if (FALSE) {
        dataset <- "sim-test-inf2"
        scen <- 1
        rep <- 1
        parents_only <- TRUE
    }

    name <- str_glue("scen-{scen}-{rep}")

    pa_dir <- str_glue("datasets/{dataset}/gfx/pred_accs")
    if (!dir.exists(pa_dir)) {
        message("- mkdir ", pa_dir)
        dir.create(pa_dir)
    }

    res <- readRDS(str_glue("datasets/{dataset}/results/{name}.rds"))
    params <- res$params

    x <- res$popn[, .SD, .SDcols = patterns("^(id|sdp|sus|inf|tol|end)")]
    if (parents_only) x <- x[sdp != "progeny"]
    setnames(x, str_replace_all(names(x),
                                c("susceptibility" = "sus",
                                  "infectivity"    = "inf",
                                  "latency"        = "lat",
                                  "detectability"  = "det",
                                  "tolerance"      = "tol",
                                  "endurance"      = "end",
                                  "BV" = "g", "GV" = "g", "EV" = "e")))
    cols <- c("id", "sdp", "sus_g", "inf_g", "tol_g", "sus_e", "inf_e", "tol_e")
    setcolorder(x, cols, skip_absent = TRUE)

    if (any(cols %notin% names(x))) {
        return(NULL)
    }

    rf <- str_glue("datasets/{dataset}/data/{name}-out/summary_inf.rds")
    if (!file.exists(rf)) {
        message(str_glue("Warning: missing '{rf}'"))
        return(c(sus = 0, inf = 0, tol = 0))
    }
    y <- readRDS(rf)$popn[, .SD, .SDcols = patterns("^(id|sdp|sus|inf|tol|end)")]
    if (parents_only) y <- y[sdp != "progeny"]

    y <- y[, map(.SD, mean), .(id, sdp), .SDcols = patterns("_[ge]$")]
    setcolorder(y, cols, skip_absent = TRUE)

    xy <- merge(melt(x[, .SD, .SDcols = -patterns("_e")], id.vars = c("id", "sdp"), value.name = "x"),
                melt(y[, .SD, .SDcols = -patterns("_e")], id.vars = c("id", "sdp"), value.name = "y"))
    # xy <- merge(melt(x, "id", value.name = "x"),
    #             melt(y, "id", value.name = "y"))

    qx <- xy[, quantile(x, 0.8) |> as.numeric()]
    qy <- xy[, quantile(y, 0.8) |> as.numeric()]

    xy[, topx := fcase(x >= quantile(x, 0.8) & y >= quantile(y, 0.8), "TP",
                       x <  quantile(x, 0.8) & y <  quantile(y, 0.8), "TN",
                       x >= quantile(x, 0.8) & y <  quantile(y, 0.8), "FN",
                       x <  quantile(x, 0.8) & y >= quantile(y, 0.8), "FP"), variable]

    p <- ggplot(xy, aes(x = x, y = y, group = variable, colour = topx)) +
        geom_point() +
        geom_abline(linetype = "dashed", linewidth = 0.5) +
        geom_hline(yintercept = 0, colour = "black") +
        geom_vline(xintercept = 0, colour = "black") +
        geom_smooth(mapping = aes(x, y),
                    method = lm, se = FALSE, fullrange = TRUE,
                    colour = "red", linewidth = 0.5) +
        scale_colour_manual("Top 20%",
                            breaks = c("TP", "TN", "FN", "FP"),
                            labels = c("True Positive", "True Negative",
                                       "False Negative", "False Positive"),
                            # values = c("green3", "blue3", "red2", "red4")) +
                            values = brewer.pal(n = 4, name = "Set1")) +
        labs(x = "True BVs",
             y = "Estimated BVs") +
        # coord_fixed() +
        facet_wrap(vars(variable),
                   ncol = 3,
                   # scales = "free",
                   labeller = labeller(
                       variable = c(sus_g = "Genetic Susceptibility",
                                    inf_g = "Genetic Infectivity",
                                    tol_g = "Genetic Endurance",
                                    sus_e = "Environmental Susceptibility",
                                    inf_e = "Environmental Infectivity",
                                    tol_e = "Environmental Endurance"))) +
        theme_bw() +
        theme(legend.position = "bottom")
    p

    title_plt <- ggplot() +
        labs(title = with(params, str_glue("Dataset: '{dataset} / s{scenario}-{replicate}'")),
             subtitle = params$description |>
                 str_split_1(", ") |>
                 str_subset("coverage|convergence", negate = TRUE) |>
                 str_flatten_comma()) +
             theme_classic() +
             theme(plot.title = element_text(size = 22),
              plot.subtitle = element_text(size = 16))

    plt <- plot_grid(title_plt, p,
                     ncol = 1, rel_heights = c(0.5, 3))

    ggsave(str_glue("{pa_dir}/{name}-pa{parents}.pdf",
                    parents = if (parents_only) "_parents" else ""),
           plt, width = 15, height = 5)
    plt
}

pred_accs <- function(dataset = "sim-base-inf", scen = 1, rep = 1, parents_only = TRUE, method = "kendall") {
    if (FALSE) {
        dataset <- "sim-base-inf"
        scen <- 1
        rep <- 1
        parents_only <- TRUE
        method <- "kendall"
    }

    name <- str_glue("scen-{scen}-{rep}")

    x <- readRDS(str_glue("datasets/{dataset}/results/{name}.rds"))$popn
    x[, names(.SD) := NULL, .SDcols = !patterns("^(id|sdp|sus|inf|tol|end)")]
    if (parents_only) x <- x[sdp != "progeny"]
    x[, sdp := NULL]

    setnames(x, str_replace_all(names(x),
                                c("susceptibility" = "sus",
                                  "infectivity"    = "inf",
                                  "latency"        = "lat",
                                  "detectability"  = "det",
                                  "tolerance"      = "tol",
                                  "endurance"      = "end",
                                  "BV" = "g", "GV" = "g", "EV" = "e")))
    cols <- c("id", "sus_g", "inf_g", "tol_g", "sus_e", "inf_e", "tol_e")
    setcolorder(x, cols, skip_absent = TRUE)

    if (any(cols %notin% names(x))) {
        return(c(sus = 0, inf = 0, tol = 0))
    }

    rf <- str_glue("datasets/{dataset}/data/{name}-out/summary_inf.rds")
    if (!file.exists(rf)) {
        message(str_glue("Warning: missing '{rf}'"))
        return(c(sus = 0, inf = 0, tol = 0))
    }
    y <- readRDS(rf)$popn[, .SD, .SDcols = patterns("^(id|sdp|sus|inf|tol|end)")]
    if (parents_only) y <- y[sdp != "progeny"]
    y[, sdp := NULL]

    y <- y[, map(.SD, mean), id, .SDcols = patterns("_[ge]$")]
    setcolorder(y, cols, skip_absent = TRUE)

    if (str_detect("top", method)) {
        pc2p <- function(x) 1 - x / 100
        method <- str_remove(method, "top") |> as.numeric() |> pc2p()
    }

    if (is.character(method)) {
        if (method %notin% c("pearson", "spearman", "kendall")) {
            method <- "kendall"
        }
        c(sus = cor.test(x$sus_g, y$sus_g, method = method)$estimate |> unname(),
          inf = cor.test(x$inf_g, y$inf_g, method = method)$estimate |> unname(),
          tol = cor.test(x$tol_g, y$tol_g, method = method)$estimate |> unname())
    } else if (is.numeric(method)) {
        p <- method |> max(0.001) |> min(0.999)

        x[, `:=`(top_sus = sus_g >= quantile(sus_g, p),
                 top_inf = inf_g >= quantile(inf_g, p),
                 top_tol = tol_g >= quantile(tol_g, p))]

        y[, `:=`(top_sus = sus_g >= quantile(sus_g, p),
                 top_inf = inf_g >= quantile(inf_g, p),
                 top_tol = tol_g >= quantile(tol_g, p))]

        c(sus = mean(x$top_sus & y$top_sus) / (1 - p),
          inf = mean(x$top_inf & y$top_inf) / (1 - p),
          tol = mean(x$top_tol & y$top_tol) / (1 - p))
    }
}

if (FALSE) {
    dataset <- "sim-test-inf2"

    sr <- list.files(str_glue("datasets/{dataset}/results")) |>
        str_sort(numeric = TRUE) |> str_remove_all("scen-|.rds") |>
        str_split("-") |> map(as.integer)
    scens <- map_int(sr, 1)
    reps  <- map_int(sr, 2)

    methods <- c("kendall", "pearson", "spearman", "top20")
    parents_only <- TRUE

    PAs <- map(methods, \(method) {
        map2(scens, reps, possibly(
            ~ pred_accs(dataset, .x, .y, parents_only, method)
        )) |>
            map(as.list) |> rbindlist(idcol = "id")
    }) |>
        map(~ {
            .x[, `:=`(scen = scens[id], rep = reps[id])] |>
                setcolorder(c("id", "scen", "rep"))
            .x[, map(.SD, ~ round(100 * mean(.x), 1)),
               .(scen = str_c("s", scen)),
               .SDcols = patterns("sus|inf|tol")]
        }) |>
        setNames(methods)

    saveRDS(PAs, str_glue("datasets/{dataset}/meta/PAs{parents}.rds",
                          parents = if (parents_only) "-parents" else ""))
}
