{
    library(data.table)
    library(stringr)
    library(purrr)
}

duplicate_trials <- function(base = "fb_1_rpw", n = 2) {
    x <- readRDS(str_glue("fb_data/{base}.rds"))

    x2 <- rbind(
        x[is.na(trial)],
        x[trial == 1][rep(seq_len(.N), n)],
        x[trial == 2][rep(seq_len(.N), n)]
    )
    x2[, id := .I]
    x2[sdp == "progeny", group := ceiling(seq(.N) / 25)]

    base2 <- str_replace(base, "fb", str_c("fb", n))
    f <- str_glue("fb_data/{base2}.rds")
    message(str_glue("- writing '{f}'"))
    saveRDS(x2, f)

    x2
}

x2 <- duplicate_trials("fb_1_rpw", 2)
