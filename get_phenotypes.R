{
    library(data.table)
    library(stringr)
    library(purrr)
    library(HDInterval)
}

get_phenotypes <- function(dataset = "fb-test", scen = 2, rep = 1) {
    if (FALSE) {
        dataset <- "fb-test"
        scen <- 2
        rep <- 1
    }

    f <- str_glue("datasets/{dataset}/data/scen-{scen}-{rep}-out/trace_combine.tsv")
    x <- fread(f)
    x <- x[, .SD, .SDcols = !patterns("^G_|state")]

    x[, `:=`(cov_G_si = r_G_si * sqrt(cov_G_ss * cov_G_ii),
             cov_G_st = r_G_st * sqrt(cov_G_ss * cov_G_tt),
             cov_G_it = r_G_it * sqrt(cov_G_ii * cov_G_tt),
             cov_E_si = r_E_si * sqrt(cov_E_ss * cov_E_ii),
             cov_E_st = r_E_st * sqrt(cov_E_ss * cov_E_tt),
             cov_E_it = r_E_it * sqrt(cov_E_ii * cov_E_tt))]

    x[, `:=`(cov_P_ss = cov_G_ss + cov_E_ss,
             cov_P_ii = cov_G_ii + cov_E_ii,
             cov_P_tt = cov_G_tt + cov_E_tt,
             cov_P_si = cov_G_si + cov_E_si,
             cov_P_st = cov_G_st + cov_E_st,
             cov_P_it = cov_G_it + cov_E_it)]

    x[, `:=`(r_P_si = cov_P_si / sqrt(cov_P_ss * cov_P_ii),
             r_P_st = cov_P_st / sqrt(cov_P_ss * cov_P_tt),
             r_P_it = cov_P_it / sqrt(cov_P_ii * cov_P_tt))]

    x[, `:=`(h2_ss = cov_G_ss / cov_P_ss,
             h2_ii = cov_G_ii / cov_P_ii,
             h2_tt = cov_G_tt / cov_P_tt)]

    x[, c("cov_G_si", "cov_G_st", "cov_G_it",
          "cov_E_si", "cov_E_st", "cov_E_it",
          "cov_P_si", "cov_P_st", "cov_P_it") := NULL]

    cols <- c(str_subset(names(x), "_G_"),
              str_subset(names(x), "_E_"),
              str_subset(names(x), "_P_"),
              str_subset(names(x), "h2_"))

    setcolorder(x, cols)

    data.table(
        parameter = names(x),
        mean   = x[, .SD |> map(mean) |> map(round, 2)] |> unlist(),
        median = x[, .SD |> map(median) |> map(round, 2)] |> unlist(),
        lower  = x[, .SD |> map(hdi) |> map("lower") |> map(round, 2)] |> unlist(),
        upper  = x[, .SD |> map(hdi) |> map("upper") |> map(round, 2)] |> unlist()
    )
}

get_phenotypes("fb-test", 8)
