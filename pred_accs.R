{
    library(data.table)
    library(stringr)
    library(purrr)
    library(ggplot2)
    library(cowplot)
}

pred_accs_plot <- function(dataset = "sim-base-inf", name = "scen-5-1", parents_only = TRUE) {
    if (FALSE) {
        dataset <- "sim-base-inf"; name <- "scen-5-1"; parents_only <- TRUE
    }
    
    res <- readRDS(str_glue("datasets/{dataset}/results/{name}.rds"))
    params <- res$params
    
    x <- res$popn
    x[, str_subset(names(x), "^(id|sdp|sus|inf|tol|end)", negate = TRUE) := NULL]
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
    
    rf <- str_glue("datasets/{dataset}/data/{name}-out/etc_inf.rds")
    if (!file.exists(rf)) {
        message(str_glue("Warning: missing '{rf}'"))
        return(c(sus = 0, inf = 0, tol = 0))
    }
    y <- readRDS(rf)$popn
    y[, str_subset(names(y), "^(id|sdp|sus|inf|tol|end)", negate = TRUE) := NULL]
    if (parents_only) y <- y[sdp != "progeny"]
    
    y <- y[, map(.SD, mean), .(id, sdp), .SDcols = str_subset(names(y), "_[ge]$")]
    setcolorder(y, cols, skip_absent = TRUE)
    
    xy <- merge(melt(x[, .(id, sdp, sus_g, inf_g, tol_g)], id.vars = c("id", "sdp"), value.name = "x"),
                melt(y[, .(id, sdp, sus_g, inf_g, tol_g)], id.vars = c("id", "sdp"), value.name = "y"))
    # xy <- merge(melt(x, "id", value.name = "x"),
    #             melt(y, "id", value.name = "y"))
    xy[, topx := fifelse(x >= quantile(x, 0.8) & y >= quantile(y, 0.8),
                         "correct", "wrong") |> factor(), variable]
    
    p <- ggplot(xy, aes(x = x, y = y, group = variable, colour = topx)) +
        geom_point() +
        geom_abline(linetype = "dashed", linewidth = 0.5) +
        geom_hline(yintercept = 0, colour = "black") +
        geom_vline(xintercept = 0, colour = "black") +
        geom_smooth(inherit.aes = FALSE,
                    mapping = aes(x, y),
                    method = lm, se = FALSE, linewidth = 0.5) +
        scale_colour_manual("Top 20%",
                            breaks = c("correct", "wrong"),
                            labels = c("Correct", "Wrong"),
                            values = c("green3", "tomato")) +
        labs(x = "True BVs",
             y = "Estimated BVs") +
        # coord_fixed() +
        facet_wrap(. ~ variable,
                   ncol = 3,
                   labeller = labeller(
                       variable = c(sus_g = "Genetic Susceptibility",
                                    inf_g = "Genetic Infectivity",
                                    tol_g = "Genetic Endurance",
                                    sus_e = "Environmental Susceptibility",
                                    inf_e = "Environmental Infectivity",
                                    tol_e = "Environmental Endurance"))) +
        theme_bw() +
        theme(legend.position = "bottom")
    
    title_plt <- ggplot() +
        labs(title = with(params, str_glue("Dataset: '{dataset} / s{scenario}-{replicate}'")),
             subtitle = str_remove(params$description, ", (coverage|convergence)")) +
        theme_classic() +
        theme(plot.title = element_text(size = 22),
              plot.subtitle = element_text(size = 16))
    
    plt <- plot_grid(title_plt, p,
                     ncol = 1, rel_heights = c(0.2, 1))
    
    ggsave(str_glue("datasets/{dataset}/gfx/pred_accs/{name}-pa{parents}.pdf",
                    parents = if (parents_only) "_parents" else ""),
           plt, width = 15, height = 5)
    plt
}

pred_accs <- function(dataset = "sim-base-inf", name = "scen-5-1", parents_only = TRUE, method = "kendall") {
    if (FALSE) {
        dataset <- "sim-base-inf"; name <- "scen-5-1"
        parents_only <- TRUE; method <- "kendall"
    }
    
    x <- readRDS(str_glue("datasets/{dataset}/results/{name}.rds"))$popn
    x[, str_subset(names(x), "^(id|sdp|sus|inf|tol|end)", negate = TRUE) := NULL]
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
    
    rf <- str_glue("datasets/{dataset}/data/{name}-out/etc_inf.rds")
    if (!file.exists(rf)) {
        message(str_glue("Warning: missing '{rf}'"))
        return(c(sus = 0, inf = 0, tol = 0))
    }
    y <- readRDS(rf)$popn
    if (parents_only) y <- y[sdp != "progeny"]
    y[, sdp := NULL]
    
    y[, str_subset(names(y), "^(id|sus|inf|tol|end)", negate = TRUE) := NULL]
    y <- y[, map(.SD, mean), id, .SDcols = str_subset(names(y), "_[ge]$")]
    setcolorder(y, cols, skip_absent = TRUE)
    
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
    scens <- 1:13
    PAs <- map(scens, \(scen) {
        reps <- 1:20
        out <- map(reps, ~ pred_accs("sim-base-inf", str_glue("scen-{scen}-{.x}"),
                                     parents_only = TRUE, method = "kendall"))
        out1 <- map(out, as.list) |> rbindlist(idcol = "rep")
        # out1[, map(.SD, mean)]
    }) |>
        rbindlist(idcol = "scen")
    
    PAs[, .(Scenario = str_c("s", scen),
            Susceptibility = round(100 * sus, 1),
            Infectivity    = round(100 * inf, 1),
            Tolerance      = round(100 * tol, 1))]
}
