{
    library(data.table)
    library(stringr)
    library(Matrix)
}

make_grm <- function(popn, use_grm) {
    switch(use_grm,
           "A" = make_A_matrix(popn) |>
               with(sparseMatrix(i = i, j = j, x = x)),
           "A_nz" = make_A_matrix(popn),
           "A_inv" = make_A_matrix(popn) |>
               with(sparseMatrix(i = i, j = j, x = x)) |>
               solve(),
           "A_nz_inv" =,
           "A_inv_nz" = make_A_matrix(popn) |>
               with(sparseMatrix(i = i, j = j, x = x)) |>
               solve() |>
               summary() |>
               as.data.table() |>
               setorder(i, j),
           # H matrix should be already constructed, just refer to NULL and
           # handle it as necessary
           NULL)
}

# Build numerator relationship "A" matrix from pedigree part of popn file
make_A_matrix <- function(popn) {
    message("Building the A-matrix ...")
    
    tmp <- popn[, .(expand.grid(i = id, j = id))]
    
    tmp[, `:=`(sire1 = popn[tmp$i, sire],
               dam1  = popn[tmp$i, dam],
               sire2 = popn[tmp$j, sire],
               dam2  = popn[tmp$j, dam])]
    
    # all related values should be 0 unless specified otherwise
    tmp[, value := 0]
    
    # full sibs
    tmp[sire1 == sire2 & dam1 == dam2, value := 1/2]
    
    # half sibs
    tmp[sire1 == sire2 & dam1 != dam2 & !is.na(dam1) & !is.na(dam2), value := 1/4]
    tmp[sire1 != sire2 & !is.na(sire1) & !is.na(sire2) & dam1 == dam2, value := 1/4]
    
    # parent/child
    tmp[i == sire2 | i == dam2 | j == sire1 | j == dam1, value := 1/2]
    
    # same individual (must be done last)
    tmp[i == j, value := 1]
    
    make_sparse <- TRUE
    if (make_sparse) {
        GRM <- tmp[value > 0, .(i, j, x = value)]
        # GRM <- with(tmp[value > 0], sparseMatrix(i = i, j = j, x = value))
    } else {
        GRM <- dcast(tmp, i ~ j, value.var = "value")
        GRM[, i := NULL]
        
        GRM <- as.matrix(GRM)
        dimnames(GRM) <- list(popn[, id], popn[, id])
    }
    
    
    GRM
}
