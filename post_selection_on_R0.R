{
    library(data.table)
    library(stringr)
    library(MASS)
    library(purrr)
    library(ggplot2)
    library(ggtext)
}

R0s <- 1:2 |>
    map(\(x) list.files("r0_selection",
                        full.names = TRUE,
                        pattern = str_glue("_{x}")) |>
            map(\(y) y |>
                    readRDS() |>
                    as.list() |>
                    as.data.table()) |>
            rbindlist() |>
            melt(measure.vars = c("none", "r0"),
                 variable.name = "selection",
                 value.name = "R0")) |>
    rbindlist(idcol = "cors")

R0s[, `:=`(cors = factor(cors),
           selection = factor(selection))]

r0_levels <- c("Heterogeneous popn", "Homogeneous popn")
# str_c("Correlation set to = ", c("+ 0.2", "  0.0", "- 0.2"))
setattr(R0s$cors, "levels", r0_levels)

# scen = "I05-no_FEs"; title_str <- bquote(I[0] == 5 ~ "no FEs")
scen = "I0_5-FEs"; title_str <- bquote(I[0] == 5 ~ "with FEs, Fishboost design")
data_name <- str_glue("r0_results/{scen}.rds")
# saveRDS(R0s, file = data_name)
# R0s <- readRDS(data_name)
# R0s <- R0s[1:2000]

mean_R0s <- R0s[, .(mean_R0 = round(mean(R0), 2),
                    median_R0 = round(median(R0), 2)),
                .(cors, selection)]

R0het <- R0s[cors == "Heterogeneous popn"]

epi_pc <- data.frame(x = c(1.4, 0.8),
                     y = factor(c("None", "R0")),
                     label = c("0%<br><i>R</i><sub>0</sub> &le; 1",
                               "21.2%<br><i>R</i><sub>0</sub> &le; 1"))

ggplot(data = R0het) +
    geom_boxplot(aes(x = R0, y = selection, fill = selection),
                 staplewidth = 0.5,
                 outliers = FALSE) +
    # geom_density(aes(x = R0, fill = selection),
    #              adjust = 2,
    #              alpha = 0.5) +
    geom_vline(xintercept = 1,
               linetype = "dashed") +
    scale_fill_discrete("Selection on",
                        labels = c("None",
                                   "<i>R</i><sub>0</sub>")) +
    geom_richtext(data = epi_pc,
                  aes(x = x, y = y, label = label)) +
    expand_limits(x = 0, y = 0) +
    guides(fill = guide_legend(reverse = TRUE)) +
    labs(x = "<i>R</i><sub>0</sub>",
         y = NULL,
         # title = bquote("Effect on" ~ R[0] ~ "from selection on Top 20% of EBVs")) +
         title = "Effect on <i>R</i><sub>0</sub> from selection on Top 20% of EBVs") +
    # theme_minimal() +
    theme_bw() +
    theme(plot.title = element_markdown(),
          axis.title.x = element_markdown(),
          legend.text = element_markdown(),
          legend.position = "inside",
          legend.position.inside = c(0.9, 0.8),
          legend.box.background = element_rect(size = 1))

ggplot(data = R0s) +
    geom_boxplot(aes(x = R0,
                     y = selection,
                     fill = selection),
                 staplewidth = 0.5,
                 outliers = FALSE) +
    scale_fill_discrete("Selection on",
                        labels = c("None", bquote(R[0]))) +
    expand_limits(x = 0,
                  y = 0) +
    guides(fill = guide_legend(reverse = TRUE)) +
    labs(x = bquote(R[0]),
         y = NULL,
         title = NULL) +
    facet_wrap(. ~ cors, ncol = 1) +
    # theme_minimal() +
    theme(legend.position.inside = c(0.8, 0.2))

ggsave(str_glue("r0_results/r0-{scen}-boxplot.png"),
       width = 9, height = 6)

ggplot(data = R0s) +
    geom_density(aes(x = R0, fill = selection),
                 adjust = 2,
                 alpha = 0.5) +
    geom_vline(data = mean_R0s,
               aes(xintercept = mean_R0, colour = selection),
               show.legend = FALSE) +
    # geom_boxplot(aes(x = R0, y = selection, fill = selection),
    #              staplewidth = 0.5,
    #              outliers = FALSE) +
    scale_fill_discrete("Selection on",
                        labels = c("None", bquote(R[0]))) +
    expand_limits(x = 0, y = 0) +
    guides(fill = guide_legend(reverse = TRUE)) +
    labs(x = bquote(R[0]),
         y = NULL,
         title = NULL) +
    facet_wrap(. ~ cors, ncol = 1) +
    # theme_minimal() +
    theme(legend.position.inside = c(0.8, 0.2))

ggsave(str_glue("r0_results/r0-{scen}-density.png"),
       width = 9, height = 6)

R0s2 <- copy(R0s)
R0s2[, Type := factor(
    fcase(cors == "Heterogeneous population" & selection == "None", "Heterogeneous\nPopulation\nNo Selection",
          cors == "Heterogeneous population" & selection == "R0",   "Heterogeneous\nPopulation\nSelection on R0",
          default =                                                   "Homogeneous\nPopulation\nNo Selection possible"))]

ggplot(data = R0s2) +
    # geom_density(aes(x = R0, fill = selection),
    #              adjust = 1,
    #              alpha = 0.5) +
    # geom_vline(data = mean_R0s,
    #            aes(xintercept = mean_R0, colour = selection),
    #            show.legend = FALSE) +
    geom_boxplot(aes(x = R0, y = Type, fill = Type),
                 staplewidth = 0.5,
                 outliers = FALSE) +
    expand_limits(x = 0, y = 0) +
    guides(fill = "none") +
    labs(x = bquote(R[0]),
         y = NULL,
         title = NULL) +
    # facet_wrap(. ~ cors, ncol = 1) +
    coord_flip() +
    theme_classic() +
    theme(legend.position.inside = c(0.8, 0.2))

ggsave(str_glue("r0_results/r0-{scen}-boxplot-alt.png"),
       width = 9, height = 6)


x <- fread("fb-final/data/scen-1-1-out/trace_combine.tsv")
x[, str_subset(names(x), "beta|period|_G_", negate = TRUE) := NULL]

x[, `:=`(
    cov_G_si = r_G_si * sqrt(cov_G_ss * cov_G_ii),
    cov_G_st = r_G_st * sqrt(cov_G_ss * cov_G_tt),
    cov_G_it = r_G_it * sqrt(cov_G_ii * cov_G_tt)
)]
pars <- as.list(x[, map(.SD, mean)])

sit <- c("sus", "inf", "tol")
Sigma <- with(pars, matrix(c(cov_G_ss, cov_G_si, cov_G_st,
                             cov_G_si, cov_G_ii, cov_G_it,
                             cov_G_st, cov_G_it, cov_G_tt),
                           3, 3,
                           dimnames = list(sit, sit)))
mu <- with(pars, beta * (latent_period + removal_period))

BVs <- mvrnorm(1e4, rep(0, 3), Sigma) |>
    as.data.table()
BVs[, r0 := exp(sus + inf + tol)]

q80 <- BVs[, unname(quantile(r0, 0.8))]
dens <- density(BVs$r0, bw = "SJ", adjust = 1, cut = 0)
dd <- with(dens, data.table(x, y))
dd[, q := fifelse(x < q80, "Remainder", "Top 80%")]

ggplot(dd) +
    geom_area(aes(x = x, y = y, fill = q)) +
    scale_fill_manual("Quantile",
                      values = c("royalblue", "tomato"),
                      labels = c("<80%", "top 80%"),
                      drop = FALSE) +
    xlim(0, 10) +
    labs(x = bquote(R[0]),
         y = "Density")

