{
    library(data.table)
    library(stringr)
    library(purrr)
}

tmp <- readRDS("datasets/sim-base-inf/results/scen-5-1.rds")
patch_pars <- tmp$params$priors[str_detect(parameter, "_[GE]_|detection"), parameter]

files <- list.files("datasets/sim-base-inf/results",
                    full.names = TRUE) |>
    str_sort(numeric = TRUE)

walk(files, \(f) {
    x <- readRDS(f)
    mf <- with(x$params, str_glue("datasets/{patch_dataset}/results/{patch_name}.rds"))
    m <- readRDS(mf)$params$priors
    
    x$params$priors[parameter %in% patch_pars,
                    true_val := m[parameter %in% patch_pars, true_val]]
    saveRDS(x, f)
})

files <- list.files("datasets/sim-base-inf/results", "-1.rds",
           full.names = TRUE) |>
    str_sort(numeric = TRUE)
foo <- data.table(file = files |>
                      str_remove_all("datasets/sim-base-inf/results/scen-|.rds"),
                  description = files |>
                      str_sort(numeric = TRUE) |>
                      map(readRDS) |>
                      map("params") |>
                      map_chr("description") |>
                      str_remove_all("FB_12_rpw, |GRM pedigree, |, convergence"),
                  patch_name = files |>
                      str_sort(numeric = TRUE) |>
                      map(readRDS) |>
                      map("params") |>
                      map_chr("patch_name") |>
                      str_remove_all("cen-|-1$"))
foo

mf <- "datasets/sim-base/results/scen-1-1.rds"
x <- readRDS(mf)
x$params$priors
x$params$priors[, true_val := fcase(
    # str_detect(parameter, "cov_._ss"), 1.0,
    # str_detect(parameter, "cov_._ii"), 0.0,
    # str_detect(parameter, "cov_._tt"), 0.5,
    # str_detect(parameter, "cov_"), 0.5,
    # str_detect(parameter, "r_._si"), 0,
    # str_detect(parameter, "r_._st"), -0.3,
    # str_detect(parameter, "r_._it"), 0.3,
    # str_detect(parameter, "r_._"), 0.2,
    default = true_val
)]
saveRDS(x, mf)
