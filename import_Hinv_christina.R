library(data.table)
library(Matrix)
library(stringr)

source("utils.R")

# Load GRM from files provided by Christina


# Load the curated FB data
fb12 <- readRDS("fb_data/fb_12.rds")
fb12 <- fb12[, .(id, sire, dam, sdp, trial, group)]

# IDs 326-350 were excluded due to tank failure
# IDs 901-1000 didn't exist

# Load the ped file and bind the fb data, we now have a lookup table
pr <- fread("Hinv/christina/ped_1775_ind_new.txt")
setnames(pr, c("new_id", "new_sire", "new_dam"))
pr <- rbind(data.table(new_id = rep(NA, 54)), pr, fill = TRUE)
pr[, c("id", "sire", "dam")] <- fb12[, .(id, sire, dam)]

# Find and fix pedigree failures
pr[V2 == 0 | V3 == 0]
pr[sire ==  9 & V2 != 0, unique(V2)]
pr[sire == 24 & V2 != 0, unique(V2)]
pr[dam  == 43 & V3 != 0, unique(V3)]
pr[dam  == 45 & V3 != 0, unique(V3)]
pr[V1 ==  221, `:=`(V2 = "T1",  V3 = "R25")]
pr[V1 == 1209, `:=`(V2 = "R14", V3 = "R31")]

my_map <- pr[, .(from = c(V1, V2, V3),
                 to = c(id, sire, dam))]
my_map <- unique(my_map, by = c("from", "to"))
my_map <- my_map[!is.na(from)]

# Load the Hinv matrix
# H <- fread("Hinv/christina/HinvOrig.txt")
H <- fread("Hinv/christina/HinvOrig_alpha1_beta0_new.txt")
setnames(H, c("i1", "j1", "x"))

# Map i,j values in Hinv
H[, i := my_map[from == i1, to], i1]
H[, j := my_map[from == j1, to], j1]
H <- H[, .(i, j, x)]

# This may be necessary if members Trial 1 Group 14 are still included
H <- H[!is.na(i) & !is.na(j)]

setorder(H, i, j)
H

Hinv12 <- with(H, sparseMatrix(i = i, j = j, x = x, symmetric = TRUE))


# Split into Trials 1 and 2
trial1_ids <- fb12[trial == 1, c(sort(unique(sire)), sort(unique(dam)), id)]
trial2_ids <- fb12[trial == 2, c(sort(unique(sire)), sort(unique(dam)), id)]

Hinv1 <- Hinv12[trial1_ids, trial1_ids]
Hinv2 <- Hinv12[trial2_ids, trial2_ids]



# Write files ----

for (ti in c(1, 2, 12)) {
    # Sparse Matrices
    M <- get(str_glue("Hinv{ti}"))
    dt <- as.data.table(summary(M))
    dt <- rbind(dt, dt[i != j, .(i = j, j = i, x)])
    setorder(dt, i, j)
    fwrite(dt, file = str_glue("fb_data/Hinv_nz{ti}.tsv"),
           sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    # write_sparse_matrix(get(str_glue("Hinv{ti}")),
    #                     file = str_glue("fb_data/Hinv_nz{ti}.tsv"))
    
    # Dense Matrices
    dt <- as.data.table(as.matrix(M))
    setnames(dt, as.character(seq_along(dt)))
    fwrite(dt, file = str_glue("fb_data/Hinv{ti}.tsv"),
           sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

