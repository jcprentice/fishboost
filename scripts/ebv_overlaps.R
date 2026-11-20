{
    library(data.table)
    library(stringr)
    library(purrr)
    library(ggplot2)
    library(ggridges)

    safe_source("figures/theme_natcom.R")
}

ebv_overlaps <- function(dataset = "fb-test",
                         scen = 2,
                         rep = 1,
                         sdp = "progeny",
                         type = "boxplot") {

    if (FALSE) {
        dataset <- "sim-test-inf1"
        scen <- 2
        rep <- 1
        sdp <- "progeny"
        type <- "boxplot"
    }

    {
        base_dir <- str_glue("datasets/{dataset}")
        data_dir <- str_glue("{base_dir}/data")
        gfx_dir  <- str_glue("{base_dir}/gfx")
        overlaps_dir <- str_glue("{gfx_dir}/overlaps")

        if (!dir.exists(overlaps_dir)) {
            message("- mkdir ", overlaps_dir)
            dir.create(overlaps_dir)
        }
    }

    sf <- str_glue("{data_dir}/scen-{scen}-{rep}-out/summary_inf.rds")
    x <- readRDS(sf)
    x <- x$popn[, .SD, .SDcols = patterns("state|id|sdp|_g$")]

    x[, r0_g := sus_g + inf_g + tol_g]

    medians <- x[, map(.SD, median), .(id, sdp), .SDcols = -"state"]

    sdp1 <- sdp
    sus_ids <- x[sdp == sdp1][order(sus_g), id[c(1, .N)]]
    inf_ids <- x[sdp == sdp1][order(inf_g), id[c(1, .N)]]
    tol_ids <- x[sdp == sdp1][order(tol_g), id[c(1, .N)]]
    r0_ids  <- x[sdp == sdp1][order(r0_g),  id[c(1, .N)]]

    summary <- rowwiseDT(
        trait=, Lower=, Upper=,
        "sus", sus_ids[[1]], sus_ids[[2]],
        "inf", inf_ids[[1]], inf_ids[[2]],
        "tol", tol_ids[[1]], tol_ids[[2]],
        "r0", r0_ids[[1]], r0_ids[[2]]
    )

    sus <- x[id %in% sus_ids, .(id = fifelse(id == sus_ids[[1]], "Lower", "Upper"),
                                variable = "sus",
                                value = sus_g)] |>
        _[order(id, value)]
    inf <- x[id %in% inf_ids, .(id = fifelse(id == inf_ids[[1]], "Lower", "Upper"),
                                variable = "inf",
                                value = inf_g)] |>
        _[order(id, value)]
    tol <- x[id %in% tol_ids, .(id = fifelse(id == tol_ids[[1]], "Lower", "Upper"),
                                variable = "end",
                                value = tol_g)] |>
        _[order(id, value)]
    r0 <- x[id %in% tol_ids, .(id = fifelse(id == tol_ids[[1]], "Lower", "Upper"),
                               variable = "r0",
                               value = r0_g)] |>
        _[order(id, value)]

    x1 <- rbind(sus, inf, tol, r0)
    lvls <- c("sus", "inf", "end", "r0")
    x1[, variable := ordered(variable, levels = lvls)]

    overlaps <- map_dbl(list(sus, inf, tol, r0), \(x) {
        d1 <- density(x[id == "Lower", value],
                      from = min(sus$value),
                      to = max(sus$value))
        d2 <- density(x[id == "Upper", value],
                      from = min(sus$value),
                      to = max(sus$value))
        sum(pmin(d1$y, d2$y)) / sum(d1$y, d2$y)
    })

    ol_str <- str_c(round(100 * overlaps, 1), "%") |>
        as.list() |>
        setNames(lvls)

    sdp1 <- switch(sdp, "sire" = "sires", "dam" = "dams", sdp)

    plt <- ggplot(x1) +
        aes(x = value, fill = id) +
        geom_vline(xintercept = 0,
                   linetype = "dashed") +
        {if (type == "boxplot") {
            geom_boxplot(staplewidth = 0.3,
                         outliers = FALSE)
        } else {
            geom_density(alpha = 0.5)
        }} +
        labs(fill = "Extreme",
             x = "EBV",
             y = NULL,
             title = str_glue("EBV posteriors of {sdp1} with lowest (red) and highest (blue) medians"),
             subtitle = str_glue("Dataset: {dataset}, scen {scen}-{rep}")) +
        facet_wrap(~ variable,
                   labeller = labeller(
                       variable = c(
                           "sus" = str_glue("Susceptibility (overlap: {ol_str$sus})"),
                           "inf" = str_glue("Infectivity (overlap: {ol_str$inf})"),
                           "end" = str_glue("Endurance (overlap: {ol_str$end})"),
                           "r0"  = str_glue("R0 (overlap: {ol_str$r0})")
                       ))) +
        theme_bw()
    plt

    gf <- str_glue("{overlaps_dir}/{dataset}-s{scen}-{rep}-overlap-{sdp}-{type}.pdf")
    ggsave(gf, plt,
           width = 17.8, height = 12, units = "cm")

    plt
}


if (FALSE) {
    ebv_overlaps("sim-test-inf1", 2, 1, "progeny", "boxplot")
    ebv_overlaps("sim-test-inf1", 2, 1, "progeny", "density")
    ebv_overlaps("sim-test-inf1", 2, 1, "sire", "boxplot")
    ebv_overlaps("sim-test-inf1", 2, 1, "sire", "density")
    ebv_overlaps("fb-test", 2, 1, "progeny", "boxplot")
    ebv_overlaps("fb-test", 2, 1, "progeny", "density")
    ebv_overlaps("fb-test", 2, 1, "sire", "boxplot")
    ebv_overlaps("fb-test", 2, 1, "sire", "density")
    ebv_overlaps("fb-test", 8, 1, "progeny", "boxplot")
    ebv_overlaps("fb-test", 8, 1, "progeny", "density")
    ebv_overlaps("fb-test", 8, 1, "sire", "boxplot")
    ebv_overlaps("fb-test", 8, 1, "sire", "density")
}
