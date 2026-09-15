{
    library(data.table)
    library(purrr)
    library(stringr)

    safe_source("rename_pars.R")
}

# Simulation-Based Calibration

sbc <- function(dataset = "sim-test-inf1", scen = 2) {
    if (FALSE) {
        dataset <- "sim-test-inf1"; scen <- 2
    }

    tc_files <- list.files(str_glue("datasets/{dataset}/data"),
                        "trace_combine.tsv",
                        full.names = TRUE,
                        recursive = TRUE) |>
        str_subset(str_glue("scen-{scen}-")) |>
        str_sort(numeric = TRUE)

    tc <- tc_files |>
        map(fread) |>
        rbindlist(idcol = "id") |>
        _[, str_subset(names(.SD), "^G_") := NULL]


    files <- list.files(str_glue("datasets/{dataset}/results"),
                        full.names = TRUE) |>
        str_subset(str_glue("scen-{scen}-")) |>
        str_sort(numeric = TRUE)

    x <- files |>
        map(readRDS) |>
        map("parameter_estimates") |>
        rbindlist(idcol = "id") |>
        _[!str_starts(parameter, "G_"),
          .(id, parameter, true_val, mean, median, sd, hdi95min, hdi95max)]

    tv <- x[id == 1 & !is.na(true_val), .(parameter, true_val)]

    inv_var <- function(a, b) {
        c1 <- (b^2 - a^2) / (2 * log(b / a))
        c2 <- ((b - a) / log(b / a))^2
        c1 - c2
    }

    {
        priors <- readRDS(files[[1]])$params$priors
        priors <- priors[parameter %in% x[, unique(parameter)]]
        priors[type == "uniform", prior_var := (val2 - val1)^2 / 12]
        priors[type == "inverse", prior_var := inv_var(val1, val2), .I]
        priors[str_starts(parameter, "cov_"), prior_var := 4 * (1 - 2 / pi)]
        priors[str_starts(parameter, "r_"), prior_var := 1 / (2 * 1.2 + 1)]
    }

    x[, Param := pretty_names(parameter)]

    x1 <- merge(x, priors[, .(parameter, prior_var)], by = "parameter", all.x = TRUE)

    x1[, rank := match(parameter, param_order)]
    setorder(x1, rank)
    x1[, rank := NULL]

    pb <- x1[!is.na(true_val),
      .(coverage = mean(hdi95min < true_val & true_val < hdi95max),
        z_score = mean((median - true_val) / sd) |> round(2),
        contraction = mean(1 - sd^2 / prior_var) |> round(2)),
      Param]

    pb
}


if (FALSE) {
    out <- sbc("sim-test-inf1", 1)
    out <- sbc("sim-test-inf1", 2)
}

if (FALSE) {
    ggplot(pb) +
        aes(x = z_score, y = contraction) +
        geom_point() +
        geom_vline(xintercept = 0,
                   linetype = "dashed") +
        facet_wrap(. ~ Param) +
        theme_bw()
}
