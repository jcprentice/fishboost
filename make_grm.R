make_grm <- function(pedigree) {
    message("Making the GRM ...")

    tmp <- pedigree[, .(expand.grid(id1 = id, id2 = id))]

    tmp[, c("sire1", "dam1") := pedigree[tmp$id1, c("sire", "dam")]]
    tmp[, c("sire2", "dam2") := pedigree[tmp$id2, c("sire", "dam")]]

    # all related values should be 0 unless specified otherwise
    tmp[, related := 0]
    # full sibs
    tmp[sire1 == sire2 & dam1 == dam2, related := 1/2]
    # half sibs
    tmp[sire1 == sire2 & dam1 != dam2 & !is.na(dam1) & !is.na(dam2), related := 1/4]
    tmp[sire1 != sire2 & !is.na(sire1) & !is.na(sire2) & dam1 == dam2, related := 1/4]
    # parent/child
    tmp[id1 == sire2 | id1 == dam2 | id2 == sire1 | id2 == dam1, related := 1/2]
    # same individual (must be done last)
    tmp[id1 == id2, related := 1]

    GRM <- dcast(tmp, id1 ~ id2, value.var = "related")
    GRM[, id1 := NULL]

    GRM <- as.matrix(GRM)
    # colnames(GRM) <- pedigree[, id]
    # rownames(GRM) <- pedigree[, id]
    dimnames(GRM) <- NULL

    GRM
}


# This is meant to be memory friendly, but it's just too slow. I might at some
# point see if it's possible to do via Rcpp.
make_grm2 <- function(pedigree) {
    N <- nrow(pedigree)

    GRM <- Matrix(0, N, N)

    for (i in 1:N) {
        for (j in 1:N) {
            # this is way too slow to actually use
        }
    }

    GRM
}
