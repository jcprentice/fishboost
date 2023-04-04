make_pedigree <- function(params) {
    message("Making pedigree ...")

    {
        setup    <- params$setup
        nsires   <- params$nsires
        ndams    <- params$ndams
        dpsire   <- params$dpsire
        ppdam    <- params$ppdam
        nparents <- params$nparents
        nprogeny <- params$nprogeny
        ntotal   <- params$ntotal
    }

    # if "fishboost", we already have a pedigree
    if (setup  == "fishboost") {
        message(" - copying FB pedigree")
        fb_data <- readRDS("fb_data/fb_data12.rds")
        ped <- fb_data[, .(id, sire, dam, sdp)]
        return(ped)
    } else if (setup == "fb1") {
        message(" - copying FB1 pedigree")
        fb_data <- readRDS("fb_data/fb_data1.rds")
        ped <- fb_data[, .(id, sire, dam, sdp)]
        return(ped)
    } else if (setup == "fb2") {
        message(" - copying FB2 pedigree")
        fb_data <- readRDS("fb_data/fb_data2.rds")
        ped <- fb_data[, .(id, sire, dam, sdp)]
        return(ped)
    }

    # create an empty data table
    ped <- data.table(id = 1:ntotal,
                      sire = NA_integer_, dam = NA_integer_,
                      sdp = factor("", levels = c("sire", "dam", "progeny")),
                      key = "id")

    sire_ids   <- seq.int(nsires)
    dam_ids    <- seq.int(nsires + 1L, nsires + ndams)
    parent_ids <- seq.int(nparents)
    prog_ids   <- seq.int(nparents + 1L, nparents + nprogeny)

    # set type
    ped[sire_ids, sdp := "sire"]
    ped[dam_ids,  sdp := "dam"]
    ped[prog_ids, sdp := "progeny"]

    # set sires and dams
    ped[prog_ids, sire := rep(sire_ids, each = dpsire * ppdam)]
    ped[prog_ids, dam := rep(dam_ids, each = ppdam)]

    return(ped)

    # tidy up
    # ped[, id := as.factor(id)]
    # ped[, sire := as.factor(sire)]
    # ped[, dam := as.factor(dam)]
    # ped[, sdp := as.factor(sdp)]

    # duplicate this but set parents' values to NA
    # ped2 <- copy(ped)
    # ped2[parent_ids, c("sire", "dam") := NA]

    # return(list(pedigree = ped, pedigree2 = ped2))
}
