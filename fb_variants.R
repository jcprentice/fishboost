{
    library(data.table)
    library(stringr)
    library(purrr)
}

fb <- readRDS("fb_data/fb_12.rds")

# Drop group 71
fb_drop71 <- fb[is.na(group) | group != 71]
fb_drop71[id %in% c(sire, dam) | sdp == "progeny"] |>
    saveRDS("fb_data/fb_12_drop71.rds")

fb_drop71[trial %in% c(NA, 1)][id %in% c(sire, dam) | sdp == "progeny"][, id := .I] |>
    saveRDS("fb_data/fb_1_drop71.rds")
fb_drop71[trial %in% c(NA, 2)][id %in% c(sire, dam) | sdp == "progeny"][, id := .I] |>
    saveRDS("fb_data/fb_2_drop71.rds")

# Drop group 71 & sire 22 (rpw = "Ricardo Pong-Wong")
fb_rpw <- fb_drop71[id != 22]
fb_rpw[, `:=`(sire = fcase(sire == 22, NA_integer_,
                           sire > 22, sire - 1L,
                           default = sire),
              dam = dam - 1L,
              id = .I)]
fb_rpw |>
    saveRDS("fb_data/fb_12_rpw.rds")
fb_rpw[trial %in% c(NA, 1)][id %in% c(sire, dam) | sdp == "progeny"][, id := .I] |>
    saveRDS("fb_data/fb_1_rpw.rds")
fb_rpw[trial %in% c(NA, 2)][id %in% c(sire, dam) | sdp == "progeny"][, id := .I] |>
    saveRDS("fb_data/fb_2_rpw.rds")
