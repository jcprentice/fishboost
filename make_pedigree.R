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
        msgs     <- params$msgs
    }

    # if "fishboost", we already have a pedigree)
    if (str_starts(setup, "fb")) {
        f <- str_glue("fb_data/{setup}.rds")
        if (file.exists(f)) {
            if (msgs) message("- copying FB pedigree: ", setup)
            pedigree <- readRDS(f)[, .(id, sire, dam, sdp, weight)]
            return(pedigree)
        } else {
            message("- no FB file, continuing without!")
        }
    }

    # create an empty data table
    pedigree <- data.table(id = seq_len(ntotal))

    sire_ids   <- seq_len(nsires)
    dam_ids    <- seq_len(ndams) + nsires
    prog_ids   <- seq_len(nprogeny) + nparents

    # set sires and dams
    pedigree[id > nparents,
             `:=`(sire = rep(sire_ids, each = dpsire * ppdam),
                  dam  = rep(dam_ids, each = ppdam))]

    # set type
    pedigree[, sdp := fcase(id <= nsires, "sire",
                            id <= nparents, "dam",
                            default = "progeny")]

    pedigree
}

