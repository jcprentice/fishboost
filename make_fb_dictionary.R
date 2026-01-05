library(data.table)
library(purrr)

make_fb_dictionary <- function() {
    
    # Import data from XLSX ----
    
    # These are the ids of the columns we want, the names they have in the
    # Masterfile, and what we're going to rename them
    cols1 <- rowwiseDT(
        id=, old=, new=,
        1,  "Original ID", "fb_id",
        # 2,  "Individual recoded", "fb_id",
        18, "Original Sire ID", "fb_sire",
        19, "Sire_rec", "sire",
        20, "Original Dam ID", "fb_dam",
        21, "Dam_rec", "dam",
        3,  "Trial", "trial",
        4,  "Box", "group"
    )
    
    # Adjust the path here if necessary
    fb <- as.data.table(
        readxl::read_excel(
            "fb_data/Fishboost Masterfile V4_donors+sex.xlsx",
            sheet = 1,
            range = readxl::cell_cols(1:40),
            .name_repair = ~ janitor::make_clean_names(.x)
        )
    )[, .SD, .SDcols = cols1$id]
    
    
    
    # Rename columns to short single word values, both for easier use, and for
    # consistency with the rest of the code:
    setnames(fb, cols1$new)
    
    
    # Fix types and missing data ----
    
    # Convert columns to integers if possible
    cols2 <- c("fb_id", "sire", "dam", "trial", "group")
    fb[, (cols2) := map(.SD, as.integer), .SDcols = cols2]
    
    # Add sire / dam / progeny as factor (although they're all progeny here):
    fb[, sdp := factor("progeny", levels = c("sire", "dam", "progeny"))]
    
    
    # Collect pedigree columns at start
    setcolorder(fb, c("fb_id", "fb_sire", "fb_dam", "sdp", "sire", "dam", "trial", "group"))
    
    
    
    # Fix missing pedigrees ----
    
    missing_pedigree <- fb[is.na(sire) | is.na(dam), fb_id]
    
    # These individuals (221 and 1209) were marked as "removed for survival" in
    # the Master file. I've given them the pedigree that makes most sense for
    # their position in the group (consistent with 5 of each sire / dam in the
    # group)
    fb[fb_id == 221,  `:=`(sire = 12, dam = 18, fb_sire = "FB-Y-1",   fb_dam = "FB-X-168")]
    fb[fb_id == 1209, `:=`(sire = 28, dam = 20, fb_sire = "FB-Y-112", fb_dam = "FB-X-189")]
    
    
    # Remove fish with no data ----
    
    # Trial 1 Group 14 was just broken, remove it completely
    fb <- fb[!(trial == 1 & group == 14)]
    
    
    # Reassign group numbers ----
    fb[trial == 1 & group > 14, group := group - 1L]
    fb[trial == 2, group := group + 35L]
    
    setorder(fb, fb_id)
    
    # Copy over the FB ID since we're about to start changing it
    fb[, id := fb_id]
    
    
    # Create the fb_dict file ----
    
    # Since sires and dams were not involved in the Fishboost experiments, they
    # lack rows in the data, so we'll add them to the top of the pedigree
    # section.
    fb_sires <- unique(fb[, .(fb_id = fb_sire, id = sire, sdp = "sire")][order(id)])
    fb_dams  <- unique(fb[, .(fb_id = fb_dam,  id = dam,  sdp =  "dam")][order(id)])
    
    # This has to be the largest ID, not the number of entries
    max_sire_id   <- max(fb_sires$id)
    max_dam_id    <- max(fb_dams$id)
    max_parent_id <- max_sire_id + max_dam_id
    
    cols4 <- c("id", "sire", "dam", "fb_id", "fb_sire", "fb_dam", "sdp", "trial", "group")
    
    fb_dict <- rbind(fb_sires, fb_dams, fb[, ..cols4],
                     fill = TRUE)
    
    setcolorder(fb_dict, cols4)
    
    
    # Sire, dam, and progeny IDs overlap, so assign them all unique ones
    fb_dict[sdp == "dam", id := id + max_sire_id]
    fb_dict[sdp == "progeny", `:=`(id  = id + max_parent_id,
                                   dam = dam + max_sire_id)]
    
    # Now use match() to shuffle the values down again, removing unused IDs
    fb_dict[, sire := match(sire, unique(id))]
    fb_dict[, dam  := match(dam,  unique(id))]
    fb_dict[, id   := match(id,   unique(id))]
    
    
    # We now how a fb_dict file suitable for either passing to SIRE, or using as
    # a base to create trait values and simulate epidemics.
    setorder(fb_dict, id)
    
    
    # Add in PrunedPed.txt
    
    pp <- fread("Hinv/finalAnalysis/PrunedPed.txt")
    
    remap <- match(fb_dict$fb_id, pp$fb_id)
    pp_id1   <- pp[remap, id1]
    pp_sire1 <- pp[remap, sire1]
    pp_dam1  <- pp[remap, dam1]
    
    fb_dict[, `:=`(pp_id = pp_id1,
                   pp_sire = pp_sire1,
                   pp_dam = pp_dam1)]
    
    setcolorder(fb_dict, c("sdp", "trial", "group"), after = ncol(fb_dict))
    
    # Save data ----
    
    # fwrite(fb_dict, file = "fb_data/fb_dict.csv")
    fwrite(fb_dict, file = "fb_data/fb_dict.tsv", sep = "\t")
    
    fb_dict
}
