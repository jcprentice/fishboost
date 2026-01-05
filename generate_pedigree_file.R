# Fishboost data ----
#
# Convert the fishboost data from XLSX to a format usable here. The final result
# is the "pedigree" data table.

# Load libraries ----

# library(readxl)
# library(janitor)
library(data.table)
library(purrr)


# Import data from XLSX ----

# These are the ids of the columns we want, the names they have in the
# Masterfile, and what we're going to rename them
cols1 <- rowwiseDT(
    id=,  old=,  new=,
    1, "Original ID", "orig_id",
    2,  "Individual recoded", "id",
    18, "Original Sire ID", "orig_sire",
    19, "Sire_rec", "sire",
    20, "Original Dam ID", "orig_dam",
    21, "Dam_rec", "dam",
    3,  "Trial", "trial",
    4,  "Box", "group",
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

setorder(fb, id)


# Fix types and missing data ----

# Convert columns to integers if possible
cols2 <- c("id", "sire", "dam", "trial", "group")
fb[, (cols2) := map(.SD, as.integer), .SDcols = cols2]

# Add sire / dam / progeny as factor (although they're all progeny here):
fb[, sdp := factor("progeny", levels = c("sire", "dam", "progeny"))]


# Collect pedigree columns at start
setcolorder(fb, c("id", "sire", "dam", "sdp"))




# Fix missing pedigrees ----

missing_pedigree <- fb[is.na(sire) | is.na(dam), orig_id]

# These individuals (221 and 1209) were marked as "removed for survival" in the
# Master file. I've given them the pedigree that makes most sense for their
# position in the group (consistent with 5 of each sire / dam in the group)
fb[orig_id == 221,  `:=`(sire = 12, dam = 18, orig_sire = "FB-Y-1", orig_dam = "FB-X-168")]
fb[orig_id == 1209, `:=`(sire = 28, dam = 20, orig_sire = "FB-Y-112", orig_dam = "FB-X-189")]

# In case we instead want to proceed with these individuals removed
# fb <- fb[!(id %in% missing_pedigree)]


# Remove fish with no data ----
#
# Fish with no data: IDs 796, 870, and 1362 have no recorded data (see Anacleto
# paper), and were discarded from analysis. Assume just missing values?

# fb <- fb[!id %in% c(796, 870, 1362)]


# Trial 1 Group 14 was just broken, remove it completely
fb <- fb[!(trial == 1 & group == 14)]


# Reassign group numbers ----
fb[trial == 1 & group > 14, group := group - 1L]
fb[trial == 2, group := group + 35L]


# Create the fb_data file ----

# Since sires and dams were not involved in the Fishboost experiments, they lack
# rows in the data, so we'll add them to the top of the pedigree section.
fb_sires <- data.table(id = sort(unique(fb$sire)), sdp = "sire")
fb_dams  <- data.table(id = sort(unique(fb$dam)),  sdp = "dam")

# This has to be the largest ID, not the number of entries
nsires   <- max(fb_sires$id)
ndams    <- max(fb_dams$id)
nparents <- nsires + ndams

pedigree <- rbind(fb_sires, fb_dams, fb, fill = TRUE)
# need to fix column order again
setcolorder(pedigree, c("id", "sire", "dam", "sdp"))


pedigree[, orig_id := as.character(orig_id)]
orig_sire_id <- pedigree[match(pedigree[sdp == "sire", id], pedigree[, sire]), orig_sire]
orig_dam_id  <- pedigree[match(pedigree[sdp == "dam",  id], pedigree[,  dam]), orig_dam]
pedigree[sdp == "sire", orig_id := orig_sire_id]
pedigree[sdp == "dam",  orig_id := orig_dam_id]


# Sire, dam, and progeny IDs overlap, so assign them all unique ones
pedigree[sdp == "dam", id := id + nsires]
pedigree[sdp == "progeny", `:=`(id  = id + nparents,
                                dam = dam + nsires)]

# Now use match() to shuffle the values down again, removing unused IDs
pedigree[, `:=`(sire = match(sire, unique(id)),
                dam  = match(dam,  unique(id)),
                id   = match(id,   unique(id)))]

setorder(pedigree, id)


# Save data ----

# Save so we don't need to source this script next time unless anything changes
# saveRDS(fb_data, file = "fb_data/fb_12.rds")
fwrite(pedigree, file = "pedigree_for_ricardo.tsv", sep = "\t", quote = TRUE)

