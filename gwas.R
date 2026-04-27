{
    library(data.table)
    library(stringr)
    library(purrr)
}

# Load and sample EBV files ----

gwas <- function(dataset = "fb-test", scen = 1, rep = 1) {
    if (FALSE) {
        dataset <- "fb-test"
        scen <- 1
        rep <- 1
    }


    f <- str_glue("datasets/{dataset}/data/scen-{scen}-{rep}-out/etc_inf.rds")
    ebvs <- readRDS(f)$popn[, .(state, id, sdp, sus_g, inf_g, tol_g, sus_e, inf_e, tol_e)]
    ebvs[, `:=`(sus_p = sus_g + sus_e,
                inf_p = inf_g + inf_e,
                tol_p = tol_g + tol_e)]

    trials <- readRDS(f)$popn[sdp == "progeny", unique(trial)]
    if (length(trials) == 1) {
        fb12 <- readRDS("fb_data/fb_12_rpw.rds")

        sires <- intersect(fb12[sdp == "sire", id], fb12[sdp == "progeny" & trial %in% trials, sire])
        all_sires <- fb12[sdp == "sire", id]
        dams <- intersect(fb12[sdp == "dam", id], fb12[sdp == "progeny" & trial %in% trials, dam])
        all_dams <- fb12[sdp == "dam", id]
        progeny <- fb12[sdp == "progeny" & trial %in% trials, id]
        all_progeny <- fb12[sdp == "progeny", id]

        nstates <- last(ebvs$state)

        ebvs[sdp == "sire",    id := rep(sires,   nstates)]
        ebvs[sdp == "dam",     id := rep(dams,    nstates)]
        ebvs[sdp == "progeny", id := rep(progeny, nstates)]

        # Missing rows
        mr <- rbind(
            data.table(sdp = "sire",    id = setdiff(all_sires,   sires)),
            data.table(sdp = "dam",     id = setdiff(all_dams,    dams)),
            data.table(sdp = "progeny", id = setdiff(all_progeny, progeny))
        )
        mr2 <- mr[rep(seq(.N), nstates)]
        mr2[, state := rep(seq_len(nstates), each = nrow(mr))]
        ebvs2 <- rbind(ebvs, mr2, fill = TRUE)

        setorder(ebvs2, state, id)
        ebvs <- ebvs2
        rm(ebvs2, mr, mr2)
    }


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


    # Create EBV files

    BVs <- expand.grid(c("g", "e", "p"), c("sus", "inf", "tol")) |>
        rev() |> apply(1, str_flatten, "_")

    out <- map(BVs, \(bv) {
        ebv <- ebvs[, mget(c("state", "id", bv))]
        x <- dcast(ebv, id ~ state, value.var = bv)[new_map]
        setnames(x, modify_if(names(x), ~ .x != "id", ~ str_c("sample", .x)))
        x[, id := ids$name]

        x[, `:=`(mean = apply(.SD, 1, mean),
                 sd = apply(.SD, 1, sd)),
          .SDcols = -1] |>
            setcolorder(c("id", "mean", "sd"))
    }) |>
        setNames(BVs)

    gwas_dir <- str_glue("../gwas/{dataset}-s{scen}-{rep}")
    if (!dir.exists(gwas_dir)) {
        message("- mkdir ", gwas_dir)
        dir.create(gwas_dir)
    }

    walk(BVs, \(bv) {
        f <- str_glue("{gwas_dir}/{x}.csv",
                      x = str_replace_all(bv, c("_g$" = "_gen",
                                                "_e$" = "_env",
                                                "_p$" = "_pt")))
        fwrite(out[[bv]], f)
    })
}
