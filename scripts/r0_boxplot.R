{
    library(data.table)
    library(purrr)
    library(ggplot2)
    library(ggtext)

    source("get_R0.R")
}

files <- list.files("datasets/sim-r0/data", "summary_sim",
                    full.names = TRUE, recursive = TRUE)

x <- map(files, readRDS)

types <- ordered(rep(c("none", "sus", "inf", "tol", "r0"), each = length(x) / 5),
                 levels = c("none", "sus", "inf", "tol", "r0"))

# Average traits
{
    popn <- map(x, "popn") |>
        map(~ .x[sdp == "progeny"]) |>
        rbindlist(idcol = "scen")
    popn <- popn[group == 1]
    popn[, `:=`(selection = types[scen],
                r0_g = sus_g + inf_g + tol_g,
                r0_e = sus_e + inf_e + tol_e)]
    popn[, r0_p := r0_g + r0_e]
    setcolorder(popn, "selection", after = "scen")
}

{
    pop2 <- popn[, .(p_inf = mean(!is.na(Tinf))),
                 .(selection, group)] |>
        dcast(... ~ selection, value.var = "p_inf")
    pop2
}

{
    foo <- popn[, map(.SD, mean), .(scen, selection), .SDcols = patterns("_[ge]")]
    foo
}

# Range of traits
{
    foo1 <- popn[, .(min_sus = min(sus_g), max_sus = max(sus_g), mean_sus = mean(sus_g),
                     min_inf = min(inf_g), max_inf = max(inf_g), mean_inf = mean(inf_g),
                     min_tol = min(tol_g), max_tol = max(tol_g), mean_tol = mean(tol_g),
                     min_r0 = min(r0_g), max_r0 = max(r0_g), mean_r0 = mean(r0_g)),
                 .(scen, selection)]
    foo1[, names(.SD) := map(.SD, round, 2), .SDcols = patterns("min|max|mean")]
    foo1
}

{
    r0s <- rbindlist(list(
        popn[, .(type = "overall", group = 1, state = 1,
                 R0 = sum(generation == 2, na.rm = TRUE) /
                     sum(generation == 1, na.rm = TRUE)), .(selection)],
        popn[, .(type = "group", state = 1,
                 R0 = sum(generation == 2, na.rm = TRUE) /
                     sum(generation == 1, na.rm = TRUE)), .(selection, group)],
        popn[, .(type = "group_rep",
                 R0 = sum(generation == 2, na.rm = TRUE) /
                     sum(generation == 1, na.rm = TRUE)), .(selection, group, state)]),
        fill = TRUE)
    r0s[, type := ordered(type, levels = c("overall", "group", "group_rep"))]
    setorder(r0s, type, selection)
    r0s
}

p1 <- ggplot(r0s,
             aes(x = selection,
                 y = R0,
                 fill = selection)) +
    geom_boxplot(width = 0.5,
                 staplewidth = 0.2,
                 outliers = FALSE,
                 show.legend = FALSE) +
    geom_hline(yintercept = 1,
               linetype = "dashed") +
    scale_y_continuous(limits = ~ range(.x, 0)) +
    scale_x_discrete(labels = levels(types)) +
    scale_fill_brewer(palette = "Set1") +
    labs(x = NULL,
         y = "R<sub>0</sub>",
         title = "R<sub>0</sub> vs selection") +
    facet_wrap(. ~ type) +
    theme_classic() +
    theme(plot.title = element_markdown(),
          axis.title.y = element_markdown())
print(p1)

ggsave("gfx/R0_selection_all.png", p1,
       width = 10, height = 6)

p2 <- ggplot(r0s[type == "group_rep"],
             aes(x = selection,
                 y = R0,
                 # group = as.character(group),
                 fill = selection)) +
    geom_boxplot(width = 0.5,
                 staplewidth = 0.2,
                 outliers = FALSE,
                 show.legend = FALSE) +
    geom_hline(yintercept = 1,
               linetype = "dashed") +
    # coord_flip() +
    scale_y_continuous(limits = ~ range(.x, 0)) +
    scale_x_discrete(labels = c(
        "none" = "No Selection",
        "sus" = "Susceptibility",
        "inf" = "Infectivity",
        "tol" = "Endurance",
        "r0" = "R0"
    )) +
    scale_fill_brewer(palette = "Set1") +
    labs(x = NULL,
         y = "R<sub>0</sub>",
         title = "R<sub>0</sub> vs selection") +
    facet_wrap(vars(group),
               labeller = labeller(
                   group = c("1" = "Trial 1",
                             "2" = "Trial 1+2"))) +
    theme_bw() +
    theme(plot.title = element_markdown(),
          axis.title.y = element_markdown())
print(p2)

ggsave("gfx/R0_selection_trials12.png", p2,
       width = 6, height = 6)
