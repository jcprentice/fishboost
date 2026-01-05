library(data.table)

get_fb_map <- function() {
    
    # Extract the fb data ----
    cols <- data.table(
        id = c(1, 2, 18, 19, 20, 21, 3, 4),
        name = c("orig_id", "id", "orig_sire", "sire", "orig_dam", "dam", "trial", "group")
    )
    
    fb <- as.data.table(
        readxl::read_excel(
            "fb_data/Fishboost Masterfile V4_donors+sex.xlsx",
            sheet = 1,
            range = readxl::cell_cols(1:40),
            .name_repair = ~ janitor::make_clean_names(.x)
        )
    )[, .SD, .SDcols = cols$id] |>
        setnames(cols$name)
    
    # Remove the two failed groups (Tr1/14, Tr2/36)
    fb <- fb[!(trial == 1 & group == 14) & !(trial == 2 & group == 36)]
    fb[, `:=`(sdp = "progeny",
              trial = NULL, group = NULL)]
    
    # Remove sire FB-Y-92
    fb[orig_sire == "FB-Y-92", `:=`(orig_sire = NA, sire = NA)]
    
    # We now have the FB dataset with just the progeny, and group 14 removed.
    # Next we want to add in the parents, they already have IDs, so don't modify
    # those for the moment.
    
    fb_sires <- data.table(id = sort(unique(fb$sire)), sdp = "sire")
    fb_dams  <- data.table(id = sort(unique(fb$dam)),  sdp = "dam")
    fb_data  <- rbind(fb_sires, fb_dams, fb, fill = TRUE)
    
    # setcolorder(fb_data, c("id", "sire", "dam", "sdp", "fb_id", "fb_sire", "fb_dam"))
    setcolorder(fb_data, c("id", "sire", "dam", "sdp", "orig_id", "orig_sire", "orig_dam"))
    
    # Now we have the FB data set with parents & progeny, but the IDs are still
    # overlapping, and they have gaps.
    
    max_sire_id <- max(fb_sires$id)
    max_dam_id  <- max(fb_dams$id)
    nparents <- max_sire_id + max_dam_id
    fb_data[sdp == "dam", id := id + max_sire_id]
    fb_data[sdp == "progeny", `:=`(id = id + nparents,
                                   dam = dam + max_sire_id)]
    
    # Now the IDs are all unique, but there are still gaps. We can remove them
    # using `match` to shuffle all the numbers down.
    
    fb_data[, sire := match(sire, unique(id))]
    fb_data[, dam  := match(dam,  unique(id))]
    fb_data[, id   := match(id,   unique(id))]
    
    # Now fill in the values for orig_id for sires and dams (currently NA)
    
    sires <- fb_data[match(fb_data[sdp == "sire", id], sire), orig_sire]
    dams  <- fb_data[match(fb_data[sdp == "dam",  id], dam),  orig_dam]
    fb_data[, orig_id := as.character(orig_id)]
    fb_data[sdp == "sire", orig_id := sires]
    fb_data[sdp == "dam",  orig_id := dams]
    
    setindex(fb_data, NULL)
    
    # Load the pruned pedigree which has the genotype names in the order they
    # appears in Hinv.txt, and the corresponding Fishboost Masterfile names.
    
    pp <- fread("Hinv_drop71/sequence/PrunedPed.txt") |>
        setnames(c("id1", "sire1", "dam1", "fb_id", "fb_sire", "fb_dam"))
    
    # We now have all the original IDs and their positions in pedrec. So match
    # the FB ID in fb_data to the FB ID in the pp file to get the map.
    
    id_map1 <- match(fb_data$orig_id, pp$fb_id)
    id_map2 <- match(pp$fb_id, fb_data$orig_id)
    
    list(fb_data = fb_data,
         pp_to_fb = id_map1,
         fb_to_pp = id_map2)
}
    
