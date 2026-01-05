# Fishboost data ----
#
# Convert the fishboost data from XLSX to a format usable here. The final result
# is the "fb_data" data table.
#
# For simulations purposes, we need "id" and "group" for the epidemic, and
# "sire" and "dam" for generating pedigree and trait values.
#
# If just using directly with SIRE 2.0, we need the "Tinf" and "Tdeath" values.
# There should be two ways of obtaining these, so we get both to check
# consistency.
#
# Note that we only know the time individuals become infectious, if there is a
# latent incubation period, then we are dealing with an SEIR model, and so we
# actually know "Tsym", not "Tinf". We'll assume these are the same unless we
# know better (e.g. "donor" individuals that show symptoms several days after
# the experiment begins should have Tinf = 0).
#
# I have been using R's data.table package for performance reasons, but this can
# be easily converted into a data frame or a tibble if preferred.


# Load libraries ----

# library(readxl)
# library(janitor)
library(data.table)
library(purrr)

import_fishboost_data <- function() {
    
    # Import data from XLSX ----
    
    # These are the ids of the columns we want, the names they have in the
    # Masterfile, and what we're going to rename them
    cols1 <- rowwiseDT(
        id=, old=, new=,
        # 1,  "Original ID", "orig_id",
        2,  "Individual recoded", "id",
        18, "Original Sire ID", "orig_sire",
        19, "Sire_rec", "sire",
        20, "Original Dam ID", "orig_dam",
        21, "Dam_rec", "dam",
        3,  "Trial", "trial",
        4,  "Box", "group",
        9,  "Weight start (g)", "weight",
        8,  "Donor (D) / Receptor (R)", "donor",
        # 11, "Description of sympthoms", "symptoms",
        12, "Day of first symptoms", "Tsym",
        13, "Day of death", "Tdeath",
        14, "Parasites", "parasites",
        35, "Days onset symptoms (susceptibility)", "sus",
        34, "Day start to end (resilience)", "res",
        36, "Days onset to death (tolerance)", "tol",
        37, "Removed (yes / no)", "removed"
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
    
    # Donor / Recipient, "D" refers to donors (Exposed / Infectious), "R" are
    # initially Susceptible.
    fb[, donor := fifelse(donor == "D", 1L, 0L)]
    
    # Parasites "YES" = true, "slaughtered" or "no" = false
    fb[, parasites := str_starts(str_to_lower(parasites), "yes")]
    
    # Removed "yes" / "no"
    fb[, removed := removed == "yes"]
    
    # Add sire / dam / progeny as factor (although they're all progeny here):
    fb[, sdp := factor("progeny", levels = c("sire", "dam", "progeny"))]
    
    
    # Collect pedigree columns at start
    setcolorder(fb, c("id", "sire", "dam", "sdp"))
    
    
    # The following columns all contain NA's, therefore are incorrectly
    # interpreted as character vectors. Need to convert to numeric.
    cols3 <- c("Tsym", "Tdeath", "sus", "res", "tol")
    fn <- function(x) as.numeric(fifelse(x == "NA", "", x))
    fb[, (cols3) := map(.SD, fn), .SDcols = cols3]
    
    
    # Fix dates ----
    
    # Dates have been interpreted as offsets from an origin date (which seems to
    # be 1899-12-30), so need to convert into R style dates.
    date0 <- as.Date("1899-12-30")
    
    # Trial no., start, and end dates
    # 1, 2014-10-02, 2015-01-14
    # 2, 2015-01-16, 2015-06-25
    trial_dates <- c(as.Date("2014-10-02"), as.Date("2015-01-16"))
    
    # Convert Tsym and Tdeath from dates to days since start of trial.
    fb[, `:=`(Tsym = as.numeric(as.Date(Tsym, date0) - trial_dates[trial]),
              Tdeath = as.numeric(as.Date(Tdeath, date0) - trial_dates[trial]))]
    
    
    ## Dates notes as errors in the Masterfile ----
    
    # These IDs are noted in the Masterfile to have had their dates of infection
    # and death switched, also resulting in "tol" < 0.
    ids_noted_in_masterfile <- c(70, 158, 368, 466, 569, 619, 658)
    
    # Test for ourselves which dates are switched:
    ids_dates_switched <- fb[Tdeath < Tsym, id]
    # check this matches:
    setequal(ids_dates_switched, ids_noted_in_masterfile)
    setdiff(ids_dates_switched, ids_noted_in_masterfile)
    # This shows us that ID 229 also needs fixing
    
    # To correct this, swap "Tsym" / "Tdeath", "res" / "sus", and negate "tol"
    fb[id %in% ids_dates_switched, `:=`(Tsym = Tdeath, Tdeath = Tsym)]
    fb[id %in% ids_dates_switched, `:=`(sus = Tsym, res = Tdeath, tol = Tdeath - Tsym)]
    
    
    # Individual 858 is reported to have started showing signs *before* the
    # start of the experiment, while 158 and 229 show signs at t = 0. We fix
    # this by assuming symptoms were shown at t = 1.
    fb[id %in% c(158, 229, 858), `:=`(Tsym = 1, sus = 1, tol = res - 1)]
    
    
    
    
    # "Tsym" and "Tdeath" should equal "sus" and "res". We need to double check
    # that is always the case, and fix it where it isn't.
    fb[, `:=`(check_Tsym = Tsym - sus,
              check_Tdeath = Tdeath - res)]
    inconsistent_ids <- fb[check_Tsym != 0 | check_Tdeath != 0, id]
    fb[id %in% inconsistent_ids]
    
    # This identifies some date inconsistencies for IDs 50, 901, 1189, and 1744.
    # Check those entries in the Masterfile, assume Tsym / Tdeath take
    # precedence unless it's clearly otherwise
    fb[id %in% inconsistent_ids, `:=`(sus = Tsym, res = Tdeath, tol = Tdeath - Tsym)]
    fb[, c("check_Tsym", "check_Tdeath") := NULL]
    
    
    ## Introduce Tinf ----
    
    # Note: we know that the Donors were infected at the start of the
    # experiment, since Tsym is (almost) never 0, that suggests there is an
    # incubation period, and Tinf occurs earlier.
    
    fb[donor == 1, Tinf := 0]
    
    
    # Fix missing pedigrees ----
    
    missing_pedigree <- fb[is.na(sire) | is.na(dam), id]
    
    # These individuals (221 and 1109) were marked as "removed for survival" in
    # the Master file. I've given them the pedigree that makes most sense for
    # their position in the group (consistent with 5 of each sire / dam in the
    # group)
    fb[id == 221,  `:=`(sire = 12, dam = 18)]
    fb[id == 1109, `:=`(sire = 28, dam = 20)]
    
    # In case we instead want to proceed with these individuals removed
    # fb <- fb[!(id %in% missing_pedigree)]
    
    
    # Remove fish with no data ----
    #
    # Fish with no data: IDs 796, 870, and 1362 have no recorded data (see
    # Anacleto paper), and were discarded from analysis. Assume just missing
    # values?
    
    # fb <- fb[!id %in% c(796, 870, 1362)]
    
    
    # Trial 1 Group 14 was just broken, remove it completely
    fb <- fb[!(trial == 1 & group == 14)]
    
    
    # Reassign group numbers ----
    fb[trial == 1 & group > 14, group := group - 1L]
    fb[trial == 2, group := group + 35L]
    
    
    # Create the fb_data file ----
    
    # Since sires and dams were not involved in the Fishboost experiments, they
    # lack rows in the data, so we'll add them to the top of the pedigree
    # section.
    fb_sires <- data.table(id = sort(unique(fb$sire)), sdp = "sire")
    fb_dams  <- data.table(id = sort(unique(fb$dam)),  sdp = "dam")
    
    # This has to be the largest ID, not the number of entries
    nsires   <- max(fb_sires$id)
    ndams    <- max(fb_dams$id)
    nparents <- nsires + ndams
    
    cols4 <- c("id", "sire", "dam", "sdp", "trial", "group", "weight", "donor",
               "Tinf", "Tsym", "Tdeath", "parasites")
    
    fb_data <- rbind(fb_sires, fb_dams, fb[, ..cols4],
                     fill = TRUE)
    # need to fix column order again
    setcolorder(fb_data, c("id", "sire", "dam", "sdp"))
    
    
    # Sire, dam, and progeny IDs overlap, so assign them all unique ones
    fb_data[sdp == "dam", id := id + nsires]
    fb_data[sdp == "progeny", `:=`(id  = id + nparents,
                                   dam = dam + nsires)]
    
    # Now use match() to shuffle the values down again, removing unused IDs
    fb_data[, `:=`(sire = match(sire, unique(id)),
                   dam  = match(dam,  unique(id)),
                   id   = match(id,   unique(id)))]
    
    # Might be handy to keep track of how parent IDs have been reassigned. I'm
    # not saving them, but they're generated here in case we need them later on.
    fb_sires$new_id <- fb_data[sdp == "sire", id]
    fb_dams$new_id  <- fb_data[sdp == "dam",  id]
    fb_parent_ids <- rbind(fb_sires, fb_dams)
    setnames(fb_parent_ids, "id", "orig_id")
    setcolorder(fb_parent_ids, c("sdp", "orig_id", "new_id"))
    
    
    
    # We now how a fb_data file suitable for either passing to SIRE, or using as
    # a base to create trait values and simulate epidemics.
    setorder(fb_data, id)
    
    # In case it's necessary later, create a map from original to new IDs
    id_maps <- fb_data[, .(sdp, id)]
    
    setnames(id_maps, "id", "new")
    
    id_maps[sdp == "sire",    original := unique(fb$sire)]
    id_maps[sdp == "dam",     original := unique(fb$dam)]
    id_maps[sdp == "progeny", original := fb$id]
    
    setcolorder(id_maps, c("sdp", "original", "new"))
    
    # Save data ----
    
    # Save so we don't need to source this script next time unless anything changes
    saveRDS(fb_data, file = "fb_data/fb_12.rds")
    # fwrite(fb_data, file = "fb_data/fb_12.tsv", sep = "\t", quote = TRUE)
    
    fb_data
}
