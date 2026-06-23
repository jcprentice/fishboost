library(data.table)
library(stringr)
library(purrr)

files <- list.files("datasets", "scen.*.rds",
                    recursive = TRUE, full.names = TRUE) |>
    str_sort(numeric = TRUE)

rename_tp <- function(x) {
    str_replace_all(x, c("latent_period"    = "LP",
                         "detection_period" = "DP",
                         "removal_period"   = "RP"))
}

walk(files, \(f) {
    x <- readRDS(f)
    flag <- FALSE

    if ("params" %in% names(x) && any(str_detect(names(x$params), "latent_period"))) {
        flag <- TRUE
        names(x$params) <- rename_tp(names(x$params))
    }

    if (any(str_detect(x$parameter_estimates$parameter, "latent_period"))) {
        flag <- TRUE
        x$parameter_estimates[, parameter := rename_tp(parameter)]
    }

    if (any(str_detect(x$params$priors$parameter, "latent_period"))) {
        flag <- TRUE
        x$params$priors[, parameter := rename_tp(parameter)]
    }

    if (flag) {
        message(f)
        saveRDS(x, f)
    }
})


tfiles <- list.files("datasets", "trace",
                    recursive = TRUE, full.names = TRUE) |>
    str_sort(numeric = TRUE)

walk(tfiles, \(tf) {
    x <- fread(tf)
    flag <- FALSE

    if (any(str_detect(names(x), "latent_period"))) {
        flag <- TRUE
        setnames(x, rename_tp)
    }

    if (flag) {
        message(tf)
        fwrite(x, tf, sep = "\t")
    }
})
