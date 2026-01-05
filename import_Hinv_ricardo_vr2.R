library(data.table)
library(stringr)
library(Matrix)

# Load GRM from files provided by Ricardo

# Import ----

# Read the matrix, but ignore the top line which simply gives the no. of ids.
Hinv <- fread("Hinv/vr2/Hinv_ricardo_vr2.txt", skip = 1)
setnames(Hinv, c("i", "j", "x"))
setorder(Hinv, "i", "j")

# Convert to (sparse) matrix
Hinv12 <- with(Hinv, sparseMatrix(i = i, j = j, x = x, symmetric = TRUE))

# Useful value
N <- max(Hinv$i)

# Load FB data
fb12 <- readRDS("fb_data/fb_12.rds")

# Clip and sort ----

# Our goal here is to obtain a map that takes the values in x$V4 (old_id) to
# 1:1829. Note: sires = 1-29, dams = 30-54, progeny = 55-1829.

# Skip lines that are just dim size and comments
pr <- fread("Hinv/vr2/pedrec.txt", skip = 3)

setnames(pr, c("new_id", "new_sire", "new_dam", "old_id", "old_sire", "old_dam"))


# Add in original IDs. Group 14 (original IDs 326-350) were removed and the
# rest were shuffled down. Also IDs for Trial 2 start at 1001, not 901.
fb12[sdp == "progeny", original_id := c(1:325, 351:900, 1001:1900)]

# clean_id is the original ID in the masterfile, ignoring sires and dams.
pr[str_detect(old_id, "T|R", negate = TRUE),clean_id := as.integer(old_id)]

# Assign to pr the sires and dams used by fb (remove 1's, as those are missing)
matches <- match(pr$clean_id, fb12$original_id)
matches[matches == 1] <- NA

pr[, c("fb_id", "fb_sire", "fb_dam")] <- fb12[matches, .(id, sire, dam)]

new_sires <- pr[new_sire != -1, sort(unique(new_sire))]
new_dams  <- pr[new_dam  != -1, sort(unique(new_dam))]
pr[new_id %in% new_sires, "fb_id"] <- pr$fb_sire[match(new_sires, pr$new_sire)]
pr[new_id %in% new_dams,  "fb_id"] <- pr$fb_dam[match(new_dams,   pr$new_dam)]

# We have a missing sire
pr[str_detect(old_id, "T|R") & is.na(fb_id)]
pr[old_sire == 0 & !is.na(clean_id)]
setdiff(1:1829, pr$fb_id)
# fb_id = 22, old_id = T47, new_id = 49

pr[49, fb_id := 22]


# Now pr$fb_id is the map we want.
id_map <- pr[!is.na(fb_id), fb_id]

id_map
# Use the ID map we obtain from fix_pedrec.R to reorganise rows and columns
# to match the expected layout. This includes removing group 14.
Hinv12 <- Hinv12[id_map, id_map]


if (FALSE) {
    # Manual method, does not reorder individuals except 1 and 2
    
    # Handle individuals 221 and 1209 who have been moved to 2 and 1
    ids <- c(3:223, 2, 224:1211, 1, 1212:N)
    Hinv12 <- Hinv12[ids, ids]

    # Clip out group 14 (ids 380:404)
    not_group14 <- setdiff(1:N, 380:404)
    Hinv12 <- Hinv12[not_group14, not_group14]
    
    # Shift all higher ids down by 25
    Hinv[i > 404, i := i - 25]
    Hinv[j > 404, j := j - 25]
}


# Split into Trials 1 and 2
fb12 <- readRDS("fb_data/fb_12.rds")
trial1_ids <- fb12[trial == 1, c(sort(unique(sire)), sort(unique(dam)), id)]
trial2_ids <- fb12[trial == 2, c(sort(unique(sire)), sort(unique(dam)), id)]

Hinv1 <- Hinv12[trial1_ids, trial1_ids]
Hinv2 <- Hinv12[trial2_ids, trial2_ids]


# Sparse Matrices ----

# Convert to data.table and write as tsv
for (ti in c(1, 2, 12)) {
    M <- get(str_glue("Hinv{ti}"))
    
    # Sparse Matrix
    dt <- as.data.table(summary(M))
    dt <- rbind(dt, dt[i != j, .(i = j, j = i, x)])
    setorder(dt, i, j)
    fwrite(dt, file = str_glue("fb_data/Hinv_nz{ti}.tsv"),
           sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

    # Dense Matrices
    dt <- as.data.table(as.matrix(M))
    setnames(dt, as.character(seq_along(dt)))
    fwrite(dt, file = str_glue("fb_data/Hinv{ti}.tsv"),
           sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}


# We have Hinv, can use to generate H
uninvert <- FALSE
if (uninvert) {
    H12 <- solve(Hinv12)
    
    # Clip out group 14
    for (k in 1:2) {
        IDs <- fb12[trial == k, c(sort(unique(sire)), sort(unique(dam)), id)]
        dt <- H12[IDs, IDs]
        write.table(dt, file = str_glue("fb_data/H{k}.tsv"),
                    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    }
    
    write.table(H, file = "fb_data/H12.tsv",
                sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}
