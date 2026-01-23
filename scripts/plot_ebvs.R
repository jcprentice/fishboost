{
    library(data.table)
    library(stringr)
    library(purrr)
    library(HDInterval)
    library(ggplot2)
    library(RColorBrewer)
    library(ggh4x)
}

plot_ebvs <- function(dataset = "fb-test", scen = 1, rep = 1) {
    if (FALSE) {
        dataset <- "fb-test"; scen <- 1; rep <- 1
    }
    
    {
        base_dir <- str_glue("datasets/{dataset}")
        data_dir <- str_glue("{base_dir}/data")
        res_dir  <- str_glue("{base_dir}/results")
        meta_dir <- str_glue("{base_dir}/meta")
        gfx_dir  <- str_glue("{base_dir}/gfx")
        ebvs_dir <- str_glue("{gfx_dir}/ebvs")
    }
    
    walk(c(meta_dir, gfx_dir, ebvs_dir), \(d) {
        if (!dir.exists(d)) {
             message(" - mkdir ", d)
             dir.create(d)
         }
    })
    
    etc <- str_glue("{data_dir}/scen-{scen}-{rep}-out/etc_inf.rds")
    if (file.exists(etc)) {
        popn <- readRDS(etc)$popn
    } else {
        return(NULL)
    }
    
    traits <- c("sus_g", "inf_g", "tol_g", "sus_e", "inf_e", "tol_e")
    cols <- c("state", "id", "sire", "dam", "sdp", "trial", "donor")
    popn <- popn[, mget(c(cols, traits))]
    
    plts <- list()
    
    st <- str_glue("{dataset} / s{scen}-{rep}")
    
    x1 <- popn[, .(trial = first(trial),
                   sus_g__lo = hdi(sus_g)[["lower"]], sus_g__hi = hdi(sus_g)[["upper"]], sus_g__med = median(sus_g),
                   inf_g__lo = hdi(inf_g)[["lower"]], inf_g__hi = hdi(inf_g)[["upper"]], inf_g__med = median(inf_g),
                   tol_g__lo = hdi(tol_g)[["lower"]], tol_g__hi = hdi(tol_g)[["upper"]], tol_g__med = median(tol_g),
                   
                   sus_e__lo = hdi(sus_e)[["lower"]], sus_e__hi = hdi(sus_e)[["upper"]], sus_e__med = median(sus_e),
                   inf_e__lo = hdi(inf_e)[["lower"]], inf_e__hi = hdi(inf_e)[["upper"]], inf_e__med = median(inf_e),
                   tol_e__lo = hdi(tol_e)[["lower"]], tol_e__hi = hdi(tol_e)[["upper"]], tol_e__med = median(tol_e)),
               id] |>
        melt(id.var = c("id", "trial"),
             measure.vars = measure(trait, val, sep = "__")) |>
        dcast(... ~ val) |>
        setcolorder(c("id", "trial", "trait", "lo", "med", "hi"))
    x1[, trait := factor(trait, levels = traits)]
    setorder(x1, id, trait)
    
    plts$ids <- x1 |>
        ggplot(aes(x = id, colour = trait)) +
        geom_ribbon(aes(ymin = lo, ymax = hi, fill = trait)) +
        # geom_errorbar(aes(ymin = lo, ymax = hi)) +
        geom_line(aes(y = med), linewidth = 1, colour = "grey", alpha = 0.5) +
        geom_hline(yintercept = 0, colour = "black") +
        labs(x = "Id",
             y = "Trait",
             title = "Trait estimates",
             subtitle = st) +
        # scale_colour_manual(values = rep(brewer.pal(4, "Set1"),
        #                                  length.out = length(levels(x2$id)))) +
        facet_wrap(. ~ trait,
                   ncol = 3,
                   labeller = labeller(
                       trait = c(sus_g = "Sus G", inf_g = "Inf G", tol_g = "Tol G",
                                 sus_e = "Sus E", inf_e = "Inf E", tol_e = "Tol E")
                   )) +
        theme_bw() +
        theme(legend.position = "none")
    plts$ids
    
    ggsave(str_glue("{ebvs_dir}/{dataset}-s{scen}-r{rep}-ebvs.png"),
           plts$ids, width = 9, height = 6)
    
    x11 <- x1[order(trait, med)][, id := seq(.N), trait][order(id, trait)]
    
    plts$sorted <- x11 |>
        ggplot(aes(x = id, colour = trait)) +
        geom_ribbon(aes(ymin = lo, ymax = hi, fill = trait)) +
        # geom_errorbar(aes(ymin = lo, ymax = hi)) +
        geom_line(aes(y = med), linewidth = 1, colour = "grey", alpha = 0.5) +
        geom_hline(yintercept = 0, colour = "black") +
        labs(x = "Sorted Position",
             y = "Trait",
             title = "Trait estimates sorted by median",
             subtitle = st) +
        # scale_colour_manual(values = rep(brewer.pal(4, "Set1"),
        #                                  length.out = length(levels(x2$id)))) +
        facet_wrap(. ~ trait,
                   ncol = 3,
                   labeller = labeller(
                       trait = c(sus_g = "Sus G", inf_g = "Inf G", tol_g = "Tol G",
                                 sus_e = "Sus E", inf_e = "Inf E", tol_e = "Tol E")
                   )) +
        theme_bw() +
        theme(legend.position = "none")
    plts$sorted
    
    ggsave(str_glue("{ebvs_dir}/{dataset}-s{scen}-r{rep}-ebvs-sorted.png"),
           plts$sorted, width = 9, height = 6)
    
    x2 <- popn[!is.na(sire),
               .(sus_g__lo = hdi(sus_g)[["lower"]], sus_g__hi = hdi(sus_g)[["upper"]], sus_g__med = median(sus_g),
                 inf_g__lo = hdi(inf_g)[["lower"]], inf_g__hi = hdi(inf_g)[["upper"]], inf_g__med = median(inf_g),
                 tol_g__lo = hdi(tol_g)[["lower"]], tol_g__hi = hdi(tol_g)[["upper"]], tol_g__med = median(tol_g),
                 sus_e__lo = hdi(sus_e)[["lower"]], sus_e__hi = hdi(sus_e)[["upper"]], sus_e__med = median(sus_e),
                 inf_e__lo = hdi(inf_e)[["lower"]], inf_e__hi = hdi(inf_e)[["upper"]], inf_e__med = median(inf_e),
                 tol_e__lo = hdi(tol_e)[["lower"]], tol_e__hi = hdi(tol_e)[["upper"]], tol_e__med = median(tol_e)),
               sire] |>
        melt(id.var = "sire", 
             measure.vars = measure(trait, val, sep = "__")) |>
        dcast(... ~ val) |>
        setcolorder(c("sire", "trait", "lo", "med", "hi"))

    x2[, trait := factor(trait, levels = traits)]
    x2 <- x2[order(trait, med)]
    x2[, sire := seq(.N), trait]

    plts$sires <- x2 |>
        ggplot(aes(x = sire, colour = trait)) +
        geom_errorbar(aes(ymin = lo, ymax = hi)) +
        geom_point(aes(y = med)) +
        geom_hline(yintercept = 0, colour = "black") +
        labs(x = "Sire",
             y = "Trait",
             title = "Sire Trait estimates",
             subtitle = st) +
        # scale_colour_manual(values = rep(brewer.pal(4, "Set1"),
        #                                  length.out = length(levels(x2$sire)))) +
        facet_wrap(. ~ trait,
                   ncol = 3,
                   labeller = labeller(
                       trait = c(sus_g = "Sus G", inf_g = "Inf G", tol_g = "Tol G",
                                 sus_e = "Sus E", inf_e = "Inf E", tol_e = "Tol E")
                   )) +
        theme_bw() +
        theme(legend.position = "none")
    plts$sires
    
    ggsave(str_glue("{ebvs_dir}/{dataset}-s{scen}-r{rep}-ebvs-sires1.png"),
           plts$sires, width = 9, height = 6)
    
    # Using stat_summary
    x3 <- popn[state %in% 1:10 & !is.na(sire), mget(c("sire", traits))] |>
        melt(id.vars = "sire", variable.name = "trait")
    x3[, trait := factor(trait, levels = traits)]
    
    hdi95 <- function(x) {
        x <- na.omit(x)
        y <- hdi(x)
        data.frame(y = median(x),
                   ymin = y[["lower"]],
                   ymax = y[["upper"]])
    }
    
    my_med <- function(x) {
        data.frame(y = median(x, na.rm = TRUE))
    }
    
    plts$hdi <- x3 |>
        ggplot(aes(x = sire, y = value, colour = trait)) +
        stat_summary(geom = "errorbar",
                     fun.data = hdi95,
                     show.legend = FALSE) +
        stat_summary(geom = "point",
                     fun.data = my_med,
                     show.legend = FALSE) +
        geom_hline(yintercept = 0, colour = "grey10") +
        labs(x = "Sire",
             y = "Trait",
             title = "Sire Traits 95% HDI",
             subtitle = st) +
        facet_wrap(. ~ trait,
                   ncol = 3,
                   labeller = labeller(
                       trait = c(sus_g = "Sus G", inf_g = "Inf G", tol_g = "Tol G",
                                 sus_e = "Sus E", inf_e = "Inf E", tol_e = "Tol E"))) +
        theme_bw()
    plts$hdi
    
    ggsave(str_glue("{ebvs_dir}/{dataset}-s{scen}-r{rep}-ebvs-sires_hdi.png"),
           plts$hdi, width = 9, height = 6)
    
    
    popn2 <- popn |> melt(id.vars = c("id", "sire", "sdp", "trial"),
                          measure.vars = traits,
                          variable.name = "trait")
    popn2[is.na(trial), trial := 0]
    popn2[, trial := factor(trial)]
    
    plts$hist <- popn2 |>
        ggplot(aes(x = value, group = trial, fill = trait)) +
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
                       trait = c(sus_g = "Sus G", inf_g = "Inf G", tol_g = "Tol G",
                                 sus_e = "Sus E", inf_e = "Inf E", tol_e = "Tol E"),
                       trial = c("0" = "sires",
                                 "1" = "Trial 1",
                                 "2" = "Trial 2")
                   )) +
        theme_bw() +
        theme(legend.position = "none")
    plts$hist
    
    ggsave(str_glue("{ebvs_dir}/{dataset}-s{scen}-r{rep}-ebvs-hist.png"),
           plts$hist, width = 9, height = 6)
    
    
    saveRDS(plts, str_glue("{meta_dir}/{dataset}-{scen}-{rep}-ebvs.rds"))
    
    invisible(plts)
}


plot_ebvs2 <- function(dataset = "sim-base-inf", scen = 1, rep = 1) {
    if (FALSE) {
        dataset <- "sim-base-inf"; scen <- 1; rep <- 1
    }
    
    {
        base_dir <- str_glue("datasets/{dataset}")
        data_dir <- str_glue("{base_dir}/data")
        res_dir  <- str_glue("{base_dir}/results")
        meta_dir <- str_glue("{base_dir}/meta")
        gfx_dir  <- str_glue("{base_dir}/gfx")
        ebvs_dir <- str_glue("{gfx_dir}/ebvs")
    }
    
    walk(c(meta_dir, gfx_dir, ebvs_dir),
         \(d) if (!dir.exists(d)) {
             message(" - mkdir ", d)
             dir.create(d)
         }
    )
    
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
    
    traits <- c("sus_g", "inf_g", "tol_g", "sus_e", "inf_e", "tol_e")
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
                   ncol = 3,
                   labeller = labeller(
                       trait = c(sus_g = "Sus G", inf_g = "Inf G", tol_g = "Tol G",
                                 sus_e = "Sus E", inf_e = "Inf E", tol_e = "Tol E")
                   )) +
        theme_bw() +
        theme(legend.position = "none")
    plts$progeny
    
    ggsave(str_glue("{ebvs_dir}/{dataset}-s{scen}-r{rep}-ebvs-progeny.png"),
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
                   ncol = ntraits,
                   labeller = labeller(
                       trait = c(sus_g = "Sus G", inf_g = "Inf G", tol_g = "Tol G",
                                 sus_e = "Sus E", inf_e = "Inf E", tol_e = "Tol E")
                   )) +
        theme_bw() +
        theme(legend.position = "none")
    plts$sires
    
    ggsave(str_glue("{ebvs_dir}/{dataset}-s{scen}-r{rep}-ebvs-sires.png"),
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
                       trait = c(sus_g = "Sus G", inf_g = "Inf G", tol_g = "Tol G",
                                 sus_e = "Sus E", inf_e = "Inf E", tol_e = "Tol E"),
                       trial = c("0" = "sires",
                                 "1" = "Trial 1",
                                 "2" = "Trial 2")
                   )) +
        theme_bw() +
        theme(legend.position = "none")
    plts$hist
    
    ggsave(str_glue("{ebvs_dir}/{dataset}-s{scen}-r{rep}-ebvs-hist.png"),
           plts$hist, width = 9, height = 6)
    
    
    saveRDS(plts, str_glue("{meta_dir}/{dataset}-{scen}-{rep}-ebvs2.rds"))
    
    invisible(plts)
}

if (FALSE) {
    p1 <- plot_ebvs("fb-test", 1, 1)
    p2 <- plot_ebvs2("sim-base-inf", 1, 1)
}
