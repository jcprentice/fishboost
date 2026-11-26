{
    library(data.table)
    library(stringr)
    library(purrr)
    library(ggplot2)
    library(cowplot)
    library(RColorBrewer)
    library(pipebind)
}

pred_accs_plot <- function(dataset = "sim-test-inf", scen = 2, rep = 1, sdp = "sire") {
    if (FALSE) {
        dataset <- "sim-test-inf1"
        scen <- 2
        rep <- 1
        sdp <- c("all", "sire", "progeny")[2]
    }

    name <- str_glue("scen-{scen}-{rep}")

    {
        base_dir <- str_glue("datasets/{dataset}")
        data_dir <- str_glue("{base_dir}/data")
        res_dir  <- str_glue("{base_dir}/results")
        gfx_dir  <- str_glue("{base_dir}/gfx")
        pa_dir   <- str_glue("{gfx_dir}/pred_accs")
    }

    if (!dir.exists(pa_dir)) {
        message("- mkdir ", pa_dir)
        dir.create(pa_dir)
    }

    results <- readRDS(str_glue("{res_dir}/{name}.rds"))
    params <- results$params

    sdp_ <- if (length(sdp) == 0 || "all" %in% sdp) {
        c("sire", "dam", "progeny")
    } else if ("parents" %in% sdp) {
        c("sire", "dam")
    } else {
        intersect(sdp, c("sire", "dam", "progeny"))
    }

    x <- results$popn[sdp %in% sdp_, .SD, .SDcols = patterns("(^id$|_g$)")]
    setnames(x, str_replace_all(names(x),
                                c("susceptibility" = "sus",
                                  "infectivity"    = "inf",
                                  "latency"        = "lat",
                                  "detectability"  = "det",
                                  "tolerance"      = "tol",
                                  "endurance"      = "end",
                                  "BV" = "g", "GV" = "g")))
    cols <- c("id", "sus_g", "inf_g", "tol_g")
    setcolorder(x, cols, skip_absent = TRUE)

    if (any(cols %notin% names(x))) {
        message("Warning: missing columns!")
        return(NULL)
    }

    rf <- str_glue("{data_dir}/{name}-out/summary_inf.rds")
    if (!file.exists(rf)) {
        message(str_glue("Warning: missing '{rf}'"))
        return(c(sus = 0, inf = 0, tol = 0))
    }
    y <- readRDS(rf)$popn[sdp %in% sdp_, .SD, .SDcols = patterns("(^id$|_g$)")]

    y <- y[, map(.SD, median), id, .SDcols = patterns("_g$")]
    setcolorder(y, cols, skip_absent = TRUE)

    xy <- merge(melt(x[, .SD, .SDcols = patterns("id|_g$")],
                     id.vars = "id", value.name = "x"),
                melt(y[, .SD, .SDcols = patterns("id|_g$")],
                     id.vars = "id", value.name = "y"))

    qx <- xy[, quantile(x, 0.8) |> as.numeric()]
    qy <- xy[, quantile(y, 0.8) |> as.numeric()]

    xy[, topx := fcase(x >= quantile(x, 0.8) & y >= quantile(y, 0.8), "TP",
                       x <  quantile(x, 0.8) & y <  quantile(y, 0.8), "TN",
                       x >= quantile(x, 0.8) & y <  quantile(y, 0.8), "FN",
                       x <  quantile(x, 0.8) & y >= quantile(y, 0.8), "FP"), variable]

    p <- ggplot(xy) +
        aes(x = x, y = y, group = variable, colour = topx) +
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
                    parents = if (any("sire", "dam") %in% sdp_) "_parents" else ""),
           plt, width = 15, height = 5)
    plt
}

pred_accs <- function(dataset = "sim-test-inf1",
                      scen = 2,
                      rep = 1,
                      sdp = "all",
                      method = "pearson") {
    if (FALSE) {
        dataset <- "sim-n-inf"
        scen <- 1
        rep <- 1
        sdp <- c("all", "sire", "progeny")[2]
        method <- "pearson"
    }

    name <- str_glue("scen-{scen}-{rep}")

    {
        base_dir <- str_glue("datasets/{dataset}")
        data_dir <- str_glue("{base_dir}/data")
        res_dir  <- str_glue("{base_dir}/results")
    }

    sdp_ <- if (length(sdp) == 0 || "all" %in% sdp) {
        c("sire", "dam", "progeny")
    } else if ("parents" %in% sdp) {
        c("sire", "dam")
    } else {
        intersect(sdp, c("sire", "dam", "progeny"))
    }

    x <- readRDS(str_glue("{res_dir}/{name}.rds"))$popn |>
        _[sdp %in% sdp_, .SD, .SDcols = patterns("(^id$|_g$)")]

    setnames(x, str_replace_all(names(x),
                                c("susceptibility" = "sus",
                                  "infectivity"    = "inf",
                                  "latency"        = "lat",
                                  "detectability"  = "det",
                                  "tolerance"      = "tol",
                                  "endurance"      = "tol",
                                  "BV" = "g", "GV" = "g")))
    cols <- c("id", "sus_g", "inf_g", "tol_g")
    setcolorder(x, cols, skip_absent = TRUE)

    x[, r0 := sus_g + inf_g + tol_g]

    if (any(cols %notin% names(x))) {
        return(c(sus = 0, inf = 0, tol = 0, r0 = 0))
    }

    rf <- str_glue("{data_dir}/{name}-out/summary_inf.rds")
    if (!file.exists(rf)) {
        message(str_glue("Warning: missing '{rf}'"))
        return(c(sus = 0, inf = 0, tol = 0, r0 = 0))
    }
    y <- readRDS(rf)$popn[sdp %in% sdp_, .SD,
                          .SDcols = patterns("(^id$|_g$)")]

    y <- y[, map(.SD, mean), id, .SDcols = patterns("_g$")]
    setcolorder(y, cols, skip_absent = TRUE)

    y[, r0 := sus_g + inf_g + tol_g]

    if (str_detect(method, "top")) {
        pc2p <- function(x) 1 - x / 100
        method <- str_remove(method, "top") |> as.numeric() |> pc2p()
    }

    if (is.character(method)) {
        if (method %notin% c("pearson", "spearman", "kendall")) {
            method <- "pearson"
        }
        c(sus = cor.test(x$sus_g, y$sus_g, method = method)$estimate |> unname(),
          inf = cor.test(x$inf_g, y$inf_g, method = method)$estimate |> unname(),
          tol = cor.test(x$tol_g, y$tol_g, method = method)$estimate |> unname(),
          r0  = cor.test(x$r0,    y$r0,    method = method)$estimate |> unname())
    } else if (is.numeric(method)) {
        p <- method |> max(0.001) |> min(0.999)

        x[, `:=`(top_sus = sus_g >= quantile(sus_g, p),
                 top_inf = inf_g >= quantile(inf_g, p),
                 top_tol = tol_g >= quantile(tol_g, p),
                 top_r0  = r0    >= quantile(r0,    p))]

        y[, `:=`(top_sus = sus_g >= quantile(sus_g, p),
                 top_inf = inf_g >= quantile(inf_g, p),
                 top_tol = tol_g >= quantile(tol_g, p),
                 top_r0  = r0    >= quantile(r0,    p))]

        c(sus = mean(x$top_sus & y$top_sus) / (1 - p),
          inf = mean(x$top_inf & y$top_inf) / (1 - p),
          tol = mean(x$top_tol & y$top_tol) / (1 - p),
          r0  = mean(x$top_r0  & y$top_r0)  / (1 - p))
    }
}

if (FALSE) {
    dataset <- "sim-test-inf1"

    {
        base_dir <- str_glue("datasets/{dataset}")
        res_dir  <- str_glue("{base_dir}/results")
        meta_dir <- str_glue("{base_dir}/meta")
    }

    sr <- list.files(res_dir) |>
        str_remove_all("scen-|.rds") |>
        str_sort(numeric = TRUE) |>
        str_split("-") |> map(as.integer)
    scens <- map_int(sr, 1)
    reps  <- map_int(sr, 2)

    methods <- c("pearson", "kendall", "spearman", "top20")[1]
    sdp <- c("sire", "dam", "progeny")[3]

    PAs <- map(methods, \(method) {
        PAs <- map2(scens, reps, possibly(
            ~ pred_accs(dataset, .x, .y, sdp, method)
        )) |>
            map(as.list) |> rbindlist(idcol = "id")
    }) |>
        map(~ {
            .x[, `:=`(scen = scens[id], rep = reps[id])] |>
                setcolorder(c("id", "scen", "rep"))
            .x[, map(.SD, ~ round(100 * mean(.x), 1)),
               .(scen = str_c("s", scen)),
               .SDcols = patterns("sus|inf|tol|r0")]
        }) |>
        setNames(methods)

    saveRDS(PAs, str_glue("{meta_dir}/PAs{parents}.rds",
                          parents = if (any("parents", "sire", "dam") %in% sdp) "-parents" else ""))
}

