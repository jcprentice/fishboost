{
    library(data.table)
    library(stringr)
    library(purrr)
}

# Load and sample EBV files ----

dataset <- "fb-test"; scen <- 9

f <- str_glue("datasets/{dataset}/data/scen-{scen}-1-out/etc_inf.rds")
ebvs <- readRDS(f)$popn[, .(state, id, sus_g, inf_g, tol_g, sus_e, inf_e, tol_e)]
ebvs[, `:=`(sus_pt = sus_g + sus_e,
            inf_pt = inf_g + inf_e,
            tol_pt = tol_g + tol_e)]



# Maps ----
source("get_fb_map.R")

# Load the PrunedPed file
pp <- fread("Hinv_drop71/sequence/PrunedPed.txt") |>
    setnames(c("id1", "sire1", "dam1", "fb_id", "fb_sire", "fb_dam"))

# get the map between PP and FB masterfile
maps <- get_fb_map()

dt <- maps$fb_data
dt[, `:=`(mf = pp[maps$pp_to_fb, fb_id],
          gt = pp[maps$pp_to_fb, id1],
          pp_to_fb = maps$pp_to_fb,
          fb_to_pp = maps$fb_to_pp)]


# awk '{print $1}' PrunedGenotype.txt > ids.txt
ids <- fread("../gwas/ids.txt") |> setnames("id")
ids[, name := str_c("id", id)]

# rows <- data.table(snp = str_c("snp", seq(0,1274747)))

# get the new map and check it against dt
new_map <- match(ids$id, dt$gt)
dt2 <- dt[new_map]

# Create EBV siles

{
    gts_sus <- dcast(ebvs[, .(state, id, sus_g)], id ~ state, value.var = "sus_g")[new_map]
    gts_inf <- dcast(ebvs[, .(state, id, inf_g)], id ~ state, value.var = "inf_g")[new_map]
    gts_tol <- dcast(ebvs[, .(state, id, tol_g)], id ~ state, value.var = "tol_g")[new_map]
    evs_sus <- dcast(ebvs[, .(state, id, sus_e)], id ~ state, value.var = "sus_e")[new_map]
    evs_inf <- dcast(ebvs[, .(state, id, inf_e)], id ~ state, value.var = "inf_e")[new_map]
    evs_tol <- dcast(ebvs[, .(state, id, tol_e)], id ~ state, value.var = "tol_e")[new_map]
    pts_sus <- dcast(ebvs[, .(state, id, sus_pt)], id ~ state, value.var = "sus_pt")[new_map]
    pts_inf <- dcast(ebvs[, .(state, id, inf_pt)], id ~ state, value.var = "inf_pt")[new_map]
    pts_tol <- dcast(ebvs[, .(state, id, tol_pt)], id ~ state, value.var = "tol_pt")[new_map]
}

{
    gts_sus[, id := ids$name]
    gts_inf[, id := ids$name]
    gts_tol[, id := ids$name]
    evs_sus[, id := ids$name]
    evs_inf[, id := ids$name]
    evs_tol[, id := ids$name]
    pts_sus[, id := ids$name]
    pts_inf[, id := ids$name]
    pts_tol[, id := ids$name]
}

{
    cols <- c("id", str_c("sample", seq_len(ncol(ebvs_sus) - 1L)))
    setnames(gts_sus, cols)
    setnames(gts_inf, cols)
    setnames(gts_tol, cols)
    setnames(evs_sus, cols)
    setnames(evs_inf, cols)
    setnames(evs_tol, cols)
    setnames(pts_sus, cols)
    setnames(pts_inf, cols)
    setnames(pts_tol, cols)
}

# Add in mean and SD
{
    gts_sus[, `:=`(mean = apply(.SD, 1, mean), sd = apply(.SD, 1, sd)), .SDcols = -1]
    gts_inf[, `:=`(mean = apply(.SD, 1, mean), sd = apply(.SD, 1, sd)), .SDcols = -1]
    gts_tol[, `:=`(mean = apply(.SD, 1, mean), sd = apply(.SD, 1, sd)), .SDcols = -1]
    evs_sus[, `:=`(mean = apply(.SD, 1, mean), sd = apply(.SD, 1, sd)), .SDcols = -1]
    evs_inf[, `:=`(mean = apply(.SD, 1, mean), sd = apply(.SD, 1, sd)), .SDcols = -1]
    evs_tol[, `:=`(mean = apply(.SD, 1, mean), sd = apply(.SD, 1, sd)), .SDcols = -1]
    pts_sus[, `:=`(mean = apply(.SD, 1, mean), sd = apply(.SD, 1, sd)), .SDcols = -1]
    pts_inf[, `:=`(mean = apply(.SD, 1, mean), sd = apply(.SD, 1, sd)), .SDcols = -1]
    pts_tol[, `:=`(mean = apply(.SD, 1, mean), sd = apply(.SD, 1, sd)), .SDcols = -1]
}

# Reorder
{
    setcolorder(gts_sus, c("id", "mean", "sd"))
    setcolorder(gts_inf, c("id", "mean", "sd"))
    setcolorder(gts_tol, c("id", "mean", "sd"))
    setcolorder(evs_sus, c("id", "mean", "sd"))
    setcolorder(evs_inf, c("id", "mean", "sd"))
    setcolorder(evs_tol, c("id", "mean", "sd"))
    setcolorder(pts_sus, c("id", "mean", "sd"))
    setcolorder(pts_inf, c("id", "mean", "sd"))
    setcolorder(pts_tol, c("id", "mean", "sd"))
}

fwrite(gts_sus, file = "../gwas/fb_gts_sus.csv")
fwrite(gts_inf, file = "../gwas/fb_gts_inf.csv")
fwrite(gts_tol, file = "../gwas/fb_gts_tol.csv")
fwrite(evs_sus, file = "../gwas/fb_evs_sus.csv")
fwrite(evs_inf, file = "../gwas/fb_evs_inf.csv")
fwrite(evs_tol, file = "../gwas/fb_evs_tol.csv")
fwrite(pts_sus, file = "../gwas/fb_pts_sus.csv")
fwrite(pts_inf, file = "../gwas/fb_pts_inf.csv")
fwrite(pts_tol, file = "../gwas/fb_pts_tol.csv")
