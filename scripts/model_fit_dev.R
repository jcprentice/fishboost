{
    library(data.table)
    library(stringr)
    library(purrr)
    library(ggplot2)
    library(ggtext)
    library(cowplot)
    # library(ggh4x)

    safe_source("utils.R")
}

# gg_colour_hue <- function(n) {
#     hues = seq(15, 375, length = n + 1)
#     hcl(h = hues, l = 65, c = 100)[1:n]
# }

model_fit_dev <- function(dataset = "fb-test", scens = 0, alt = "",
                          opts = c("extremes")[0]) {
    if (FALSE) {
        dataset <- "fb-test"
        dataset <- "sim-test-inf"
        scens <- 0; opts = character(0); alt = ""
    }

    message(str_glue("Calculating model fit deviation for '{dataset}'"))

    suffix <- str_c("",
                    if ("drop_outliers" %in% opts) "drop",
                    if ("extremes" %in% opts) "es",
                    sep = "-")

    if (str_1st(alt) %notin% c("", "-")) alt <- str_c("-", alt)

    km_data <- readRDS(str_glue("datasets/{dataset}/meta/km_data_ps.rds"))

    if (all(scens != 0)) {
        scens1 <- km_data |> names() |> str_remove("s") |> as.integer() |> intersect(scens)
        km_data <- km_data[scens1]
    }

    params <- map(km_data, "params")
    labels <- map_chr(params, "label")
    scenarios <- map_int(params, "scenario")

    # Set limits on AUC length
    tmax <- if (str_detect(dataset, "fb")) {
        c(t1 = 104, t2 = 160)
    } else {
        c(t1 = 200, t2 = 200)
    }

    # Extract RMS deviance for each id, family, trial / Tsign, RP
    fit <- map(km_data, \(x) {
        if (FALSE) {
            i <- 2; x <- copy(km_data[[i]])
        }

        # Set full-sib family combination
        # (split families that appeared in both trials)
        x$data[, family := .GRP, .(sire, dam, trial)]

        x1 <- x$data[, .(id, family, trial, donor, Tsign, Tdeath, RP, src)] |>
            setorder(id, family, trial)

        # If Tsign is missing, sub in Tdeath if available
        x1[src == "fb" & is.na(Tsign) & !is.na(Tdeath), Tsign := Tdeath]
        x1[, Tdeath := NULL]

        if ("extremes" %in% opts) {
            foo <- x1[id == last(id), .(family, trial, donor, Tsign)]
            # Filter out any sires with fewer than 5 non-NA values
            foo[, p := .N - sum(is.na(Tsign)), .(family, trial)]
            foo <- foo[p > 5]
            foo[is.na(Tsign), Tsign := tmax[trial]]
            foo1 <- foo[, .(mu = mean(Tsign, na.rm = TRUE)), .(family, trial)] |>
                setorder(trial, family, mu)
            donor_sires <- foo[donor == 1, sort(unique(family)), trial] |>
                split(by = "trial") |> map("V1")
            foo1[, donor := fifelse(family %in% unlist(donor_sires), 1L, 0L)]
            ids <- foo1[, .(family = family[c(which.min(mu), which.max(mu))]),
                        .(trial, donor)]
            x1 <- merge(x1, ids, by = c("family", "trial", "donor")) |>
                setcolorder("donor", after = "trial")
        }

        # This only affects simulations
        x1[, `:=`(Tsign = ceiling(Tsign), RP = ceiling(RP))]

        sv_curve <- function(x) c(0, sort(x, na.last = TRUE))

        x2 <- x1[, .(survival = seq(1, 0, length = .N + 1),
                     Tsign = sv_curve(Tsign),
                     RP    = sv_curve(RP),
                     src   = first(src)),
                 .(id, family, trial)] |>
            melt(measure.vars = c("Tsign", "RP"),
                 value.name = "time") |>
            setorder(id, family, variable) |>
            setcolorder(c("src", "id", "family", "trial", "variable", "time", "survival"))

        if (FALSE && "keep_small_groups" %notin% opts) {
            # If times are only 0 or NA, then we drop them
            x2[src == "fb", keep := any(time %notin% c(0, NA)),
               .(id, family, trial, variable)]
            families <- x2[src == "fb" & keep == TRUE, unique(family)]
            pc <- x2[src == "fb", mean(keep)]
            if (pc < 1) {
                message("- Dropping families with no infection")
                x2 <- x2[family %in% families]
            }
            x2[, keep := NULL]
        }

        # Weight points
        if (FALSE) {
            # Weight by no. of points and by 1 / time covered
            x2[, `:=`(wt1 = .N,
                      wt2 = min(max(time, na.rm = TRUE), tmax[trial])),
               .(id, family, trial, variable)]
            x2[, wt := {
                x <- (wt1 * mean(wt2)) / (mean(wt1) * wt2)
                x / mean(x)
            }]
            x2[, c("wt1", "wt2") := NULL]
        } else {
            # Weight by no. of points only
            x2[, wt := .N, .(id, family, trial, variable)]
            x2[, wt := wt / mean(wt)]
        }

        # Resample x2 at all times 0:tmax[trial]
        # Don't do anything with NA values, will nafill them later
        x3 <- x2[, {
            times <- seq(0, tmax[first(trial)])
            list(time = times,
                 survival = approx(time, survival, times,
                                   method = "constant",
                                   ties = list("ordered", max))$y |>
                     # Do we want this for all values or just sim?
                     nafill("locf"),
                 src = first(src),
                 wt = first(wt))
        },
        .(id, family, trial, variable)]

        # Extend all missing sim values
        # x3[src == "sim", survival := nafill(survival, "locf")]

        # Remove times where FB is NA
        # x3[, tmp := seq(.N), id]
        # ids_to_keep <- x3[id == last(id) & !is.na(survival), tmp]
        # x3 <- x3[tmp %in% ids_to_keep]
        # x3[, tmp := NULL]

        if (FALSE) {
            # Mean and SE across id
            x4 <- x3[, .(survival_mean = mean(survival, na.rm = TRUE),
                         survival_sem = sem(survival),
                         wt = first(wt)),
                     .(family, trial, variable, time, src)]

            x5 <- x4 |>
                dcast(... ~ src, value.var = c("survival_mean",
                                               "survival_sem",
                                               "wt"))
            x5[, survival_sem_fb := NULL]
            x5[, `:=`(dev = abs(survival_mean_sim - survival_mean_fb) * wt_fb,
                      sem = survival_sem_sim * wt_fb)]

            # Reduce to summary statistics
            x6 <- x5[, .(mad_all   = sum(dev),
                         mad_Tsign = sum(dev[variable == "Tsign"]),
                         mad_RP    = sum(dev[variable == "RP"]),

                         rms_all   = sqrt(sum(dev^2)),
                         rms_Tsign = sqrt(sum(dev[variable == "Tsign"]^2)),
                         rms_RP    = sqrt(sum(dev[variable == "RP"]^2)),

                         sem_all   = sum(sem),
                         sem_Tsign = sum(sem[variable == "Tsign"]),
                         sem_RP    = sum(sem[variable == "RP"]))] |> #,
                # .(family, trial)] |>
                melt(measure.vars = measure(type, variable,
                                            pattern = "(mad|sem|rms)_(.*)"))
            x6
        } else {
            # Use bootstrapping to obtain mean and SD

            # Split by id
            s1 <- split(x3, by = "id")
            n <- length(s1)

            x8 <- map(1:2000, \(i) {
                # Sample with replacement 1..n-1, n is the actual value
                # The first one gives the true mean
                ids <- if (i == 1) 1:n else c(sample(n-1, replace = TRUE), n)
                s3 <- rbindlist(s1[ids])

                # Mean and SE across resampled id
                x4 <- s3[, .(survival_mean = mean(survival, na.rm = TRUE),
                             wt = first(wt)),
                         .(family, trial, variable, time, src)]

                x5 <- x4 |>
                    dcast(... ~ src, value.var = c("survival_mean", "wt"))
                x5[, dev := abs(survival_mean_sim - survival_mean_fb) * wt_fb]

                # Reduce to summary statistics
                x6 <- x5[, .(mad_all   = sum(dev),
                             mad_Tsign = sum(dev[variable == "Tsign"]),
                             mad_RP    = sum(dev[variable == "RP"]),

                             rms_all   = sqrt(sum(dev^2)),
                             rms_Tsign = sqrt(sum(dev[variable == "Tsign"]^2)),
                             rms_RP    = sqrt(sum(dev[variable == "RP"]^2)))] |> #,
                    # .(family, trial)] |>
                    melt(measure.vars = measure(type, variable,
                                                pattern = "(mad|rms)_(.*)"))
                x6
            }) |>
                rbindlist(idcol = "id")
            x8
            # x8[, .(value = mean(value), sem = sem(value)), .(type, variable)]
        }
    }) |>
        rbindlist(idcol = "scen") # fit <- ...

    fit[, scen := ordered(labels[scen],
                          levels = str_sort(labels, numeric = TRUE))]
    setkey(fit, NULL)

    # fit[, `:=`(scen = factor(str_c("s", scen), levels = labels),
    #            id = format(id),
    #            family = format(family),
    #            trial = format(trial))]
    setorder(fit, scen, type, variable)

    # Save to meta dir (give that dir something to do!)
    meta_dir <- str_glue("datasets/{dataset}/meta")
    if (!dir.exists(meta_dir)) dir.create(meta_dir)
    saveRDS(fit, str_glue("{meta_dir}/fit.rds"))

    fit_summary <- fit[, .(mean = mean(value),
                           sd = sd(value)),
                       .(scen, type, variable)] |>
        dcast(scen ~ ..., value.var = c("mean", "sd"))
    fit_summary[, map_if(.SD, is.numeric, round, 1)]

    fit_wide <- fit[, .(mean = mean(value),
                        lo = hdi(value)[["lower"]],
                        hi = hdi(value)[["upper"]]),
                    .(scen, type, variable)] |>
        dcast(... ~ type, value.var = c("mean", "lo", "hi"))
    fit_wide[, map_if(.SD, is.numeric, round, 2)]

    fit2 <- copy(fit)
    fit2[, trial := fcase(str_detect(scen, "[123]"), "Trial 1",
                          str_detect(scen, "[456]"), "Trial 2",
                          str_detect(scen, "[789]"), "Trials 1+2")]
    fit2[, model := fcase(str_detect(scen, "[147]"), "SS End",
                          str_detect(scen, "[258]"), "MS End",
                          str_detect(scen, "[369]"), "Homogeneous")]
    fit2[, model := ordered(model, levels = c("SS End", "MS End", "Homogeneous"))]
    fit2 <- fit2[variable == "all" & !str_detect(scen, "[456]")]
    for_paper <- fit2[type == "mad",
                      .(mean = mean(value),
                        lo = hdi(value)[["lower"]],
                        hi = hdi(value)[["upper"]]),
                      .(scen, variable, trial, model)]

    p1 <- ggplot(for_paper,
           aes(x = model,
               y = mean,
               colour = variable)) +
        geom_point(show.legend = FALSE) +
        geom_errorbar(aes(ymin = lo, ymax = hi),
                      width = 0.2,
                      show.legend = FALSE) +
        # scale_y_continuous(limits = ~ range(.x, 0)) +
        scale_colour_discrete("Measurements",
                            breaks = c("all", "Tsign", "RP"),
                            # values = gg_colour_hue(3),
                            labels = c("Combined", "Tsign", "RP")) +
        labs(colour = "Measure",
             x = "Model",
             y = "F",
             title = "Fitness Metric F by Trial and Model with 95% HDI") +
        facet_wrap(vars(trial),
                   scales = "free") +
        theme_bw() +
        theme(axis.text = element_text(size = 8))
    p1
    ggsave("gfx/model_fit.pdf", p1,
           width = 15, height = 7.5, units = "cm")
    ggsave("gfx/model_fit.png", p1,
           width = 15, height = 7.5, units = "cm", dpi = "print")

    ggplot(fit_wide,
           aes(colour = variable)) +
        geom_point(aes(x = model, y = rms)) +
        # geom_segment(aes(x = model, xend = model,
        #                  y = rms - sem, yend = rms + sem)) +
        scale_y_continuous(limits = ~ range(.x, 0)) +
        scale_colour_discrete("Measurements",
                            breaks = c("all", "Tsign", "RP"),
                            # values = gg_colour_hue(3),
                            labels = c("Combined", "Tsign", "RP")) +
        labs(colour = "Measure",
             x = "Model",
             y = "Sum RMS Deviation",
             title = "Sum RMS Deviation by Trial and Model") +
        facet_wrap(vars(trial),
                   scales = "free_x") +
        theme_bw()

    fit_plt_mad <- ggplot(fit_wide, aes(x = scen, y = rms, colour = variable)) +
        geom_point(aes(x = scen, y = mad)) +
        geom_segment(aes(x = scen, xend = scen,
                         y = mad - sem, yend = mad + sem)) +
        scale_y_continuous(limits = ~ range(.x, 0)) +
        labs(x = "Scenario",
             y = "RMS of deviation",
             title = "RMS of deviation between data and model output",
             subtitle = str_glue("Dataset: '{dataset}'")) +
        theme_bw() +
        theme(plot.title = element_text(size = 16),
              plot.background = element_rect(fill = "white"))
    fit_plt_mad

    fit_plt_rms <- ggplot(fit_wide, aes(x = scen, y = rms, colour = variable)) +
        geom_point(aes(x = scen, y = rms)) +
        geom_segment(aes(x = scen, xend = scen,
                         y = rms - sem, yend = rms + sem)) +
        scale_y_continuous(limits = ~ range(.x, 0)) +
        labs(x = "Scenario",
             y = "RMS of deviation",
             title = "RMS of deviation between data and model output",
             subtitle = str_glue("Dataset: '{dataset}'")) +
        theme_bw() +
        theme(plot.title = element_text(size = 16),
              plot.background = element_rect(fill = "white"))
    fit_plt_rms

    mf_dir <- str_glue("datasets/{dataset}/gfx/model_fit")
    if (!dir.exists(mf_dir)) {
        message("- mkdir ", mf_dir)
        dir.create(mf_dir)
    }

    ggsave(str_glue("{mf_dir}/{dataset}-model_fit-rms_dev{suffix}{alt}.png"),
           fit_plt_rms, width = 8, height = 5)
    ggsave(str_glue("{mf_dir}/{dataset}-model_fit-mad_dev{suffix}{alt}.png"),
           fit_plt_mad, width = 8, height = 5)

    # families_plt <- ggplot(fit, aes(x = family, y = value, colour = variable)) +
    #     stat_summary(geom = "errorbar", fun.data = mean_se,
    #                  position = position_dodge(0.1),
    #                  width = 0.2) +
    #     stat_summary(geom = "point", fun = mean,
    #                  position = position_dodge(0.1)) +
    #     # geom_point(show.legend = FALSE,
    #     #            position = position_dodge(0.25)) +
    #     scale_colour_discrete("Values",
    #                         breaks = c("all", "Tsign", "RP"),
    #                         # values = gg_colour_hue(3),
    #                         labels = c("Combined", "Tsign", "RP")) +
    #     labs(x = "Family",
    #          y = "RMS of deviation",
    #          title = "RMS of deviation between data and model output, per family",
    #          subtitle = str_glue("Dataset: '{dataset}'")) +
    #     facet_wrap(vars(scen), ncol = 4) +
    #     theme_bw() +
    #     theme(legend.position = "bottom")
    # families_plt
    #
    # plt_str <- str_glue("{mf_base}_by_family{alt}.png")
    # ggsave(plt_str, families_plt, width = 12, height = 8)

    list(fit = fit_wide,
         fit_plt_mad = fit_plt_mad,
         fit_plt_rms = fit_plt_rms)
}

if (FALSE) {
    out_fb_final <- model_fit_dev("fb-final")
    out_fb_final <- model_fit_dev("fb-final2")
    # out_fb_final_drop <- model_fit_dev("fb-final")
    out_fb_lp <- model_fit_dev("fb-lp")
    out_fb_lp2 <- model_fit_dev("fb-lp2")
    # out_fb_lp_drop <- model_fit_dev("fb-lp")
    out_fb_simple <- model_fit_dev("fb-simple")
    out_fb_donors <- model_fit_dev("fb-donors")
    # out_fb_donors_drop <- model_fit_dev("fb-donors")
    out_sbi <- model_fit_dev("sim-base-inf", 0, "")
    out_fb_test <- model_fit_dev("fb-test", 0, "")
    out_fb_qtest <- model_fit_dev("fb-qtest", 0, "")
}


if (FALSE) {
    scens <- 7:9
    fit_plt <- fit[scen %in% str_c("s", scens) & type == "mad" & variable == "all"] |>
        ggplot() +
        aes(x = value, fill = scen) +
        geom_density(alpha = 0.5, position = "identity") +
        labs(fill = "Scenario",
             x = "Goodness of Fit, F",
             y = "Density") +
        theme_bw()
    fit_plt
    ggsave(str_glue("gfx/fit_{str_flatten(scens)}.pdf"), fit_plt,
           width = 8, height = 5)
}
