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
        if (msgs) message(" - copying FB pedigree: ", setup)
        fb_data <- readRDS(str_glue("fb_data/{setup}.rds"))
        ped <- fb_data[, .(id, sire, dam, sdp, weight)]
        return(ped)
    }
    
    # create an empty data table
    pedigree <- data.table(id = 1:ntotal,
                           sire = NA_integer_, dam = NA_integer_,
                           sdp = factor("", levels = c("sire", "dam", "progeny")))
    
    sire_ids   <- seq_len(nsires)
    dam_ids    <- seq(nsires + 1L, nsires + ndams)
    parent_ids <- seq_len(nparents)
    prog_ids   <- seq(nparents + 1L, nparents + nprogeny)
    
    # set type
    pedigree[sire_ids, sdp := "sire"]
    pedigree[dam_ids,  sdp := "dam"]
    pedigree[prog_ids, sdp := "progeny"]
    
    # set sires and dams
    pedigree[prog_ids, sire := rep(sire_ids, each = dpsire * ppdam)]
    pedigree[prog_ids, dam := rep(dam_ids, each = ppdam)]
    
    pedigree
}
