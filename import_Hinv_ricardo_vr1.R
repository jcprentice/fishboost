library(data.table)
library(Matrix)
library(stringr)
# library(gtools)

source("get_fb_map.R")

# Load GRM from files provided by Ricardo

# Import ----

# Read the matrix, but ignore the top line which simply gives the no. of ids.
Hinv <- fread("Hinv/Hinv.txt", skip = 1)
setnames(Hinv, c("i", "j", "x"))
setorder(Hinv, "i", "j")

# Convert to (sparse) matrix
Hinv12 <- with(Hinv, sparseMatrix(i = i, j = j, x = x, symmetric = TRUE))

# This should be 1829
N <- max(Hinv$i)
stopifnot(N == 1829)

# Our goal here is to map 1-1829 in Hinv to 1-1829 in fb_data.
# Note: sires = 1-29, dams = 30-54, progeny = 55-1829.

# Skip lines that are just dim size and comments
# pr <- fread("Hinv/pedrec.txt")
# pp <- fread("Hinv/PrunedPed.txt")

id_map <- get_fb_map()$pp_to_fb

# Check we've got everything
stopifnot(!any(is.na(id_map)))
stopifnot(all(sort(id_map) == 1:1829))


# Use the ID map we obtain to reorganise rows and columns to match the expected
# layout.
Hinv12 <- Hinv12[id_map, id_map]


# Split into Trials 1 and 2
fb12 <- readRDS("fb_data/fb_12.rds")
trial1_ids <- fb12[trial == 1, c(sort(unique(sire)), sort(unique(dam)), id)]
trial2_ids <- fb12[trial == 2, c(sort(unique(sire)), sort(unique(dam)), id)]

# Invert -> filter -> invert
Hinv1 <- solve(solve(Hinv12)[trial1_ids, trial1_ids])
Hinv2 <- solve(solve(Hinv12)[trial2_ids, trial2_ids])


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
