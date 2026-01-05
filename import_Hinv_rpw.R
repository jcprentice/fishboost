{
    library(data.table)
    library(stringr)
    library(purrr)
    library(Matrix)
}

# Maps ----

# Encoded 1-1803, gives equivalent T48,T17,R18,...,98,99
pedrec <- fread("Hinv_drop71/sequence/pedrec.txt", skip = 1)

# Encoded T48,T17,R18,...,98,99, with Fishboost IDs in V4
prunedped <- fread("Hinv_drop71/sequence/PrunedPed.txt")

fb_dict <- fread("fb_data/fb_dict.tsv")

# Open FB, add in Fishboost names, remove sire 22, shift down all IDs
fb <- readRDS("fb_data/fb_12.rds")
fb[, fb_id := fb_dict$fb_id]
fb_rpw <- fb[id != 22 & (group != 71 | is.na(group))]
fb_rpw[, `:=`(sire = fcase(sire == 22, NA_integer_,
                           sire > 22, sire - 1L,
                           default = sire),
              dam = dam - 1L,
              id = .I)]

fb_rpw[, pp_id := prunedped[, V1[match(fb_id, V4)]]]
fb_rpw[, H_id := pedrec[, V1[match(pp_id, V4)]]]
setcolorder(fb_rpw, c("id", "fb_id", "pp_id", "H_id"))

if (FALSE) {
    fbx <- fb_rpw[, .SD, .SDcols = !c("fb_id", "pp_id", "H_id")]
    fbx |>
        saveRDS("fb_data/fb_12_rpw.rds")
    fbx[trial %in% c(NA, 1)][id %in% c(sire, dam) | sdp == "progeny"][, id := .I] |>
        saveRDS("fb_data/fb_1_rpw.rds")
    fbx[trial %in% c(NA, 2)][id %in% c(sire, dam) | sdp == "progeny"][, id := .I] |>
        saveRDS("fb_data/fb_2_rpw.rds")
    rm(fbx)
}


# Sparse Matrices ----

# The data I want, encoded 1-1803
HS_inv <- fread("Hinv_drop71/sequence/Hinv.txt", skip = 1) |>
    setNames(c("i", "j", "x"))
HG_inv <- fread("Hinv_drop71/imputed/Hinv.txt", skip = 1) |>
    setNames(c("i", "j", "x"))

# Be *really* careful about the match here.
# Hinv i = 1:5 corresponds to FB id = 22,1,4,2,9
HS_inv[, `:=`(i1 = match(i, fb_rpw$H_id),
              j1 = match(j, fb_rpw$H_id))]
HS_inv[, `:=`(i = pmin(i1, j1) - 1L,
              j = pmax(i1, j1) - 1L)]
HS_inv[, c("i1", "j1") := NULL]
setorder(HS_inv, i, j)
fwrite(HS_inv, "fb_data/HS_inv_nz_12_rpw.tsv", sep = "\t")

HG_inv[, `:=`(i1 = match(i, fb_rpw$H_id),
              j1 = match(j, fb_rpw$H_id))]
HG_inv[, `:=`(i = pmin(i1, j1) - 1L,
              j = pmax(i1, j1) - 1L)]
HG_inv[, c("i1", "j1") := NULL]
setorder(HG_inv, i, j)
fwrite(HG_inv, "fb_data/HG_inv_nz_12_rpw.tsv", sep = "\t")


# Full Matrices ----

n <- last(HS_inv$i) + 1L

with(HS_inv, sparseMatrix(i = i, j = j, x = x,
                          symmetric = TRUE, index1 = FALSE)) |>
    as.matrix() |> as.data.table() |>
    set_names(1:n) |>
    fwrite(file = "fb_data/HS_inv_12_rpw.tsv", sep = "\t")

with(HG_inv, sparseMatrix(i = i, j = j, x = x,
                          symmetric = TRUE, index1 = FALSE)) |>
    as.matrix() |> as.data.table() |>
    set_names(1:n) |>
    fwrite(file = "fb_data/HG_inv_12_rpw.tsv", sep = "\t")

fb <- readRDS("fb_data/fb_12_rpw.rds")
t1_sires <- fb[trial == 1, sire |> unique() |> sort()]
t2_sires <- fb[trial == 2, sire |> unique() |> sort()]
t1_ids <- fb[sire %in% t1_sires | id %in% t1_sires, id]
t2_ids <- fb[sire %in% t2_sires | id %in% t2_sires, id]

HG <- fread("fb_data/HG_inv_12_rpw.tsv", header = TRUE) |>
    as.matrix() |> solve()

HG_inv_1 <- HG[t1_ids, t1_ids] |> solve() |> as.data.table()
fwrite(HG_inv_1, file = "fb_data/HG_inv_1_rpw.tsv")

HG_inv_2 <- HG[t2_ids, t2_ids] |> solve() |> as.data.table()
fwrite(HG_inv_2, file = "fb_data/HG_inv_2_rpw.tsv")

