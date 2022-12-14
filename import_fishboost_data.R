# Fishboost data ----
#
# Convert the fishboost data from XLSX to a format usable here. The final result
# is the "fb_traits" data table.
#
# For simulations purposes, we need "id" and "group" for the epidemic, and
# "sire" and "dam" for generating pedigree and trait values.
#
# If just using directly with SIRE 2.0, we need the "Tinf" and "Trec" values.
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


# Import data from XLSX ----

# These are the ids of the columns we want, the names they have in the
# Masterfile, and what we're going to rename them
cols <- as.data.table(dplyr::tribble(
    ~id, ~old, ~new,
    2,  "Individual recoded", "id",
    19, "Sire_rec", "sire",
    21, "Dam_rec", "dam",
    3,  "Trial", "trial",
    4,  "Box", "group",
    9,  "Weight start (g)", "weight",
    8,  "Donor (D) / Receptor (R)", "donor",
    12, "Day of first symptoms", "Tsym",
    13, "Day of death", "Trec",
    35, "Days onset symptoms (susceptibility)", "sus",
    34, "Day start to end (resilience)", "res",
    36,  "Days onset to death (tolerance)", "tol"
))

# Adjust the path here if necessary
fb <- as.data.table(
    readxl::read_excel(
        "fb_data/Fishboost Masterfile V4_donors+sex.xlsx",
        sheet = 1,
        range = readxl::cell_cols(1:40),
        .name_repair = ~ janitor::make_clean_names(.x)
    )
)[, .SD, .SDcols = cols$id]



# Rename columns to short single word values, both for easier use, and for
# consistency with the rest of the code:
setnames(fb, cols$new)

# Setting the key to "id" keeps the data sorted in a useful way
setkey(fb, id)


# Fix types and missing data ----

# Convert columns to integers if possible
cols <- c("id", "sire", "dam", "trial", "group")
fb[, (cols) := lapply(.SD, as.integer), .SDcols = cols]

# Donor / Recipient, "D" refers to donors (Exposed / Infectious), "R" are
# initially Susceptible.
fb[, donor := fifelse(donor == "D", 1L, 0L)]

# Add sire / dam / progeny as factor (although they're all progeny here):
fb[, sdp := factor("progeny", levels = c("sire", "dam", "progeny"))]


# Collect pedigree columns at start
setcolorder(fb, c(1:3, ncol(fb)))


# The following columns all contain NA's, therefore are incorrectly interpreted
# as character vectors. Need to convert to numeric.
cols <- c("Tsym", "Trec", "sus", "res", "tol")
fn <- function(x) as.numeric(fifelse(x == "NA", "", x))
fb[, (cols) := lapply(.SD, fn), .SDcols = cols]


# Fix dates ----

# Dates have been interpreted as offsets from an origin date (which seems to be
# 1899-12-30), so need to convert into R style dates.
date0 <- as.Date("1899-12-30")

# Trial no., start, and end dates
# 1, 2014-10-02, 2015-01-14
# 2, 2015-01-16, 2015-06-25
trial_dates <- c(as.Date("2014-10-02"), as.Date("2015-01-16"))

# Convert Tsym and Trec from dates to days since start of trial.
fb[, `:=`(Tsym = as.numeric(as.Date(Tsym, date0) - trial_dates[trial]),
          Trec = as.numeric(as.Date(Trec, date0) - trial_dates[trial]))]


## Dates notes as errors in the Masterfile ----

# These IDs are noted in the Masterfile to have had their dates of infection and
# death switched, also resulting in "tol" < 0.
ids_noted_in_masterfile <- c(70, 158, 368, 466, 569, 619, 658)

# Test for ourselves which dates are switched:
ids_dates_switched <- fb[Trec < Tsym, id]
# check this matches:
setequal(ids_dates_switched, ids_noted_in_masterfile)
setdiff(ids_dates_switched, ids_noted_in_masterfile)
# This shows us that ID 229 also needs fixing

# To correct this, swap "Tsym" <-> "Trec", "res" <-> "sus", and negate "tol".
fb[id %in% ids_dates_switched, `:=`(Tsym = Trec, Trec = Tsym)]
fb[id %in% ids_dates_switched, `:=`(sus = Tsym, res = Trec, tol = Trec - Tsym)]


# Individual 858 is reported to have started showing signs *before* the start of
# the experiment, while 158 and 229 show signs at t = 0. We fix this by assuming
# symptoms were shown at t = 1.
fb[id %in% c(158, 229, 858), `:=`(Tsym = 1, sus = 1, tol = res - 1)]




# "Tsym" and "Trec" should equal "sus" and "res". We need to double check that
# is always the case, and fix it where it isn't.
fb[, `:=`(check_Tsym = Tsym - sus,
          check_Trec = Trec - res)]
inconsistent_ids <- fb[check_Tsym != 0 | check_Trec != 0, id]
fb[id %in% inconsistent_ids]

# This identifies some date inconsistencies for IDs 50, 901, 1189, and 1744.
# Check those entries in the Masterfile, assume Tsym/Trec take precedence unless
# it's clearly otherwise
fb[id %in% inconsistent_ids, `:=`(sus = Tsym, res = Trec, tol = Trec - Tsym)]
fb[, c("check_Tsym", "check_Trec") := NULL]


## Introduce Tinf ----

# Note: we know that the Donors were infected at the start of the experiment,
# since "Tsym" is (almost) never 0, that suggests there is an incubation period,
# and Tinf occurs earlier.

fb[donor == 1, Tinf := 0]


# Fix missing pedigrees ----

broken_ids <- fb[is.na(sire) | is.na(dam), id]

# These individuals (221 and 1109) were marked as "removed for survival" in the
# Master file. I've given them the pedigree that makes most sense for their
# position in the group (consistent with 5 of each sire / dam in the group)
fb[id == 221,  `:=`(sire = 12, dam = 18)]
fb[id == 1109, `:=`(sire = 28, dam = 20)]

# In case we instead want to proceed with these individuals removed
# fb <- fb[!(id %in% broken_ids)]



# Create the fb_traits file ----

# Since sires and dams were not involved in the Fishboost experiments, they lack
# rows in the data, so we'll add them to the top of the pedigree section.
fb_sires <- data.table(id = sort(unique(fb$sire)), sdp = "sire")
fb_dams  <- data.table(id = sort(unique(fb$dam)),  sdp = "dam")

nsires   <- max(fb_sires$id)
ndams    <- max(fb_dams$id)
nparents <- nsires + ndams

fb_traits <- rbind(fb_sires, fb_dams,
                   fb[, .(id, sire, dam, sdp, trial, group, weight, donor, Tinf, Tsym, Trec)],
                   fill = TRUE)
# need to fix column order again
setcolorder(fb_traits, c(1,3,4,2))


# Sire, dam, and progeny IDs overlap, so assign them all unique ones
fb_traits[sdp == "dam", id := id + nsires]
fb_traits[sdp == "progeny", `:=`(id  = id + nparents,
                                 dam = dam + nsires)]

# Now use match() to shuffle the values down again, removing unused IDs
fb_traits[, `:=`(sire = match(sire, unique(id)),
                 dam  = match(dam,  unique(id)),
                 id   = match(id,   unique(id)))]

# Might be handy to keep track of how parent IDs have been reassigned. I'm not
# saving them, but they're generated here in case we need them later on.
fb_sires$new_id <- fb_traits[sdp == "sire", id]
fb_dams$new_id  <- fb_traits[sdp == "dam",  id]
fb_parent_ids <- rbind(fb_sires, fb_dams)
setnames(fb_parent_ids, "id", "orig_id")
setcolorder(fb_parent_ids, c(2, 1, 3))


# The groups are repeated, we need them to also be unique.
trial1_ngroups <- fb[trial == 1, length(unique(group))]
fb_traits[trial == 2, group := group + trial1_ngroups]


# Correct for Latent Period ----
#
# If we assume an SIR model, then Tsym is the same as Tinf, in which case donors
# need to take the lower of the two values. If we assume an SEIR model, then we
# still need a Tinf value for SIRE, so take the average Tsym - Tinf for donors
# (since this is our best estimate for the latent period), and apply that to the
# receptors

correct_for_LP <- FALSE

if (correct_for_LP) {
    # we may want to do this per trial, or only using Trial 1
    mean_lp <- fb_traits[donor == 1, mean(Tsym - Tinf, na.rm = TRUE)]
    # floor(mean_lp) to ensure rounding to nearst day
    fb_traits[donor == 0, Tinf := pmax(Tsym - floor(mean_lp), 0.0)]
}

# We now how a fb_traits file suitable for either passing to SIRE, or using as a
# base to create trait values and simulate epidemics.
setkey(fb_traits, id)



# In case it's necessary later, create a map from original IDs to the new IDs.
id_maps <- fb_traits[, .(sdp, id)]

setnames(id_maps, "id", "new")

id_maps[sdp == "sire",    original := unique(fb$sire)]
id_maps[sdp == "dam",     original := unique(fb$dam)]
id_maps[sdp == "progeny", original := fb$id]

setcolorder(id_maps, c("sdp", "original", "new"))

# Save data ----

# Save so we don't need to source this script next time unless anything changes
saveRDS(fb_traits, file = "fb_data/fb_traits.rds")
# fwrite(fb_traits, file = "fb_data/fb_traits.tsv", sep = "\t", quote = TRUE)
