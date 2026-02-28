plot_ebvs2 <- function(dataset = "sim-base-inf", scen = 1, rep = 1) {
    if (FALSE) {
        dataset <- "sim-test-inf"; scen <- 1; rep <- 1
    }

    {
        base_dir <- str_glue("datasets/{dataset}")
        data_dir <- str_glue("{base_dir}/data")
        res_dir  <- str_glue("{base_dir}/results")
        meta_dir <- str_glue("{base_dir}/meta")
        gfx_dir  <- str_glue("{base_dir}/gfx")
        ebvs_dir <- str_glue("{gfx_dir}/ebvs")
    }

    c(meta_dir, gfx_dir, ebvs_dir) |>
        discard(dir.exists) |>
        walk(~ message("- mkdir ", .x)) |>
        walk(dir.create)

    etc <- str_glue("{data_dir}/scen-{scen}-{rep}-out/etc_inf.rds")
    if (file.exists(etc)) {
        popn <- readRDS(etc)$popn
    } else {
        return(NULL)
    }

    res <- str_glue("{res_dir}/scen-{scen}-{rep}.rds")
    if (file.exists(res)) {
        truevals <- readRDS(res)$popn
    } else {
        return(NULL)
    }

    traits <- c("sus_g", "inf_g", "tol_g")
    traits <- intersect(traits, intersect(names(popn), names(truevals)))
    ntraits <- traits |> str_detect("_g") |> sum()
    cols <- c("id", "sire", "dam", "sdp", "trial", "donor")
    popn <- popn[, mget(c("state", cols, traits))]
    tv1 <- truevals[, mget(c(cols, traits))] |>
        melt(id.vars = c("id", "sire", "sdp"),
             measure.vars = traits,
             variable.name = "trait")

    plts <- list()
    st <- str_glue("{dataset} / s{scen}-{rep}")

    popn[, sdp := levels(sdp)[as.integer(sdp)]]

    x1 <- popn[, c(sire = first(sire),
                   sdp = first(sdp),
                   map(.SD, quantile, 0.025) |> setNames(str_c(traits, "__lo")),
                   map(.SD, quantile, 0.5)   |> setNames(str_c(traits, "__med")),
                   map(.SD, quantile, 0.975) |> setNames(str_c(traits, "__hi"))),
               id, .SDcols = traits] |>
        melt(id.vars = c("id", "sire", "sdp"),
             measure.vars = measure(trait, val, sep = "__")) |>
        dcast(... ~ val, value.var = "value") |>
        setcolorder(c("lo", "med", "hi"), after = "trait")
    x1[, trait := factor(trait, levels = traits)]
    setorder(x1, id, trait)

    X <- merge(tv1, x1, by = c("id", "sire", "sdp", "trait")) |>
        setnames("value", "true")

    plts$progeny <- X[sdp == "progeny"] |>
        ggplot(aes(x = true, y = med, colour = trait)) +
        geom_point(alpha = 0.5) +
        geom_hline(yintercept = 0, colour = "black") +
        geom_vline(xintercept = 0, colour = "black") +
        geom_abline(colour = "black", linetype = "dashed") +
        geom_smooth(colour = "red", linewidth = 0.5, method = lm, se = FALSE) +
        # coord_equal() +
        labs(x = "True",
             y = "Median",
             title = "Progeny TBVs vs EBVs",
             subtitle = st) +
        # scale_colour_manual(values = rep(brewer.pal(4, "Set1"),
        #                                  length.out = length(levels(x2$id)))) +
        facet_wrap(. ~ trait,
                   ncol = 2,
                   labeller = labeller(
                       trait = c(sus_g = "Sus G",
                                 inf_g = "Inf G",
                                 tol_g = "Tol G")
                   )) +
        theme_bw() +
        theme(legend.position = "none")
    plts$progeny

    ggsave(str_glue("{ebvs_dir}/{dataset}-s{scen}-{rep}-ebvs-progeny.png"),
           plts$progeny, width = 9, height = 6)


    popn2 <- merge(tv1,
                  popn |> melt(id.vars = c("id", "sire", "sdp", "trial"),
                               measure.vars = traits,
                               variable.name = "trait"),
                  by = c("id", "sire", "sdp", "trait")) |>
        setnames(c("value.x", "value.y"), c("true", "estd")) |>
        setcolorder("trial", after = "sdp")


    plts$sires <- X[sdp == "sire"] |>
        ggplot(aes(x = true, y = med, colour = trait)) +
        geom_point() +
        geom_errorbar(aes(x = true, ymin = lo, ymax = hi),
                      linewidth = 0.2) +
        geom_hline(yintercept = 0, colour = "black") +
        geom_vline(xintercept = 0, colour = "black") +
        geom_abline(colour = "black", linetype = "dashed") +
        geom_smooth(data = popn2,
                    mapping = aes(x = true, y = estd),
                    method = "lm", se = FALSE,
                    colour = "red", linewidth = 0.5) +
        # geom_smooth(method = lm, se = FALSE,
        #             colour = "green", linewidth = 0.5V) +
        # coord_equal() +
        labs(x = "True",
             y = "Median",
             title = "Sire TBVs vs EBVs",
             subtitle = st) +
        # scale_colour_manual(values = rep(brewer.pal(4, "Set1"),
        #                                  length.out = length(levels(x2$id)))) +
        facet_wrap(. ~ trait,
                   ncol = 2,
                   labeller = labeller(
                       trait = c(sus_g = "Sus G",
                                 inf_g = "Inf G",
                                 tol_g = "Tol G")
                   )) +
        theme_bw() +
        theme(legend.position = "none")
    plts$sires

    ggsave(str_glue("{ebvs_dir}/{dataset}-s{scen}-{rep}-ebvs-sires.png"),
           plts$sires, width = 9, height = 6)

    popn2[is.na(trial), trial := 0]
    popn2[, trial := factor(trial)]

    plts$hist <- popn2 |>
        ggplot(aes(x = true, group = trial, fill = trait)) +
        geom_histogram(aes(y = ..density..),
                       bins = 50) +
        stat_theodensity() +
        labs(x = "EBV",
             y = "Density",
             title = "EBV density with fitted Normal distribution",
             subtitle = st) +
        facet_grid(cols = vars(trait),
                   rows = vars(trial),
                   scales = "free_y",
                   labeller = labeller(
                       trait = c(sus_g = "Sus G",
                                 inf_g = "Inf G",
                                 tol_g = "Tol G"),
                       trial = c("0" = "sires",
                                 "1" = "Trial 1",
                                 "2" = "Trial 2")
                   )) +
        theme_bw() +
        theme(legend.position = "none")
    plts$hist

    ggsave(str_glue("{ebvs_dir}/{dataset}-s{scen}-{rep}-ebvs-hist.png"),
           plts$hist, width = 9, height = 6)


    saveRDS(plts, str_glue("{meta_dir}/{dataset}-{scen}-{rep}-ebvs2.rds"))

    invisible(plts)
}

if (FALSE) {
    p2 <- plot_ebvs2("sim-base-inf", 1, 1)
}
