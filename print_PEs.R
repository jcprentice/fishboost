{
    library(data.table)
    library(stringr)
    library(purrr)
    library(DT)

    source("rename_pars.R")
}

print_PEs <- function(dataset = "fb-test", scen = 2) {
    if (FALSE) {
        dataset <- "fb-test"
        scen <- 2
    }
    x <- readRDS(str_glue("datasets/{dataset}/results/scen-{scen}-1.rds"))

    pes <- x$parameter_estimates[!str_starts(parameter, "G"),
                                 .(parameter, mean, median, hdi95min, hdi95max)]
    p2 <- merge(pes, x$params$priors[, .(parameter, type, val1, val2)], all.x = TRUE)
    setcolorder(p2, c("type", "val1", "val2"), after = "parameter")

    p2 <- p2[match(pes$parameter, p2$parameter)]

    p2[, Parameter := pretty_names(parameter)]
    p2[, Parameter := str_replace_all(Parameter,
                                      c(" Gen" = "G",
                                        " Env" = "E",
                                        " PT"  = "P"))]

    p2[, Prior := str_glue("{x}({v1}, {v2})",
                           x = fcase(str_detect(parameter, "cov"), "Half Normal",
                                     str_detect(parameter, "^r_"), "LKJ",
                                     type == "inverse", "Jeffreys",
                                     type == "uniform", "Uniform"),
                           v1 = val1, v2 = val2)]
    p2[str_detect(Prior, "LKJ"), Prior := "LKJ(1.2)"]

    p2[, d := fcase(
        str_detect(parameter, "_[GEP]_|h2_"), 2,
        str_detect(parameter, "period"), 1,
        str_detect(parameter, "weight"), 1,
        default = 1
    )]

    p2[, `:=`(
        Mean   = round(mean, 2),  # sprintf(str_glue("%5 .{d}f", d = d), mean),
        Median = round(median, 2) # sprintf(str_glue("%5 .{d}f", d = d), median)
    )]

    # p2[, `95% HDI` := sprintf(str_glue("(%5 .{d}f, %5 .{d}f)", d = d),
    #                       hdi95min, hdi95max), .I]
    p2[, `95% HDI` := str_glue("({x}, {y})",
                               x = round(hdi95min, 2),
                               y = round(hdi95max, 2)), .I]

    p2[str_detect(parameter, "_P_|h2"), `:=`(Prior = "", `95% HDI` = "")]

    p3 <- p2[, .(Parameter, Prior, Mean, Median, `95% HDI`)]
    fwrite(p3, file = "parameter_estimates.tsv", sep = "\t")
    p3
}

print_PEs()
