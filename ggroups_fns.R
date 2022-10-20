buildA <- function(ped) {
    colnames(ped) <- c("id", "sire", "dam")
    A <- diag(nrow(ped))
    for (i in 2:nrow(ped)) {
        si <- ped[i, 2:3]
        di <- si[[2]]
        si <- si[[1]]
        if (si > 0 | di > 0) {
            if (si > 0)
                si <- which(ped$id == si)
            if (di > 0)
                di <- which(ped$id == di)
            for (j in 1:i) {
                if (j == i) {
                    if (si > 0 & di > 0)
                        A[i, j] <- 1 + A[si, di]/2
                }
                else {
                    if (si > 0)
                        A[i, j] <- A[j, i] <- A[j, si]/2
                    if (di > 0)
                        A[i, j] <- A[j, i] <- A[j, i] + A[j, di]/2
                }
            }
        }
    }
    colnames(A) <- row.names(A) <- ped$id
    return(A)
}

tabA <- function(ped) {
    curr.set <- ped[ped$sire == 0 & ped$dam == 0, ]$id
    tbA <- data.frame(id1 = curr.set, id2 = curr.set, a = 1)
    ped <- ped[!ped$id %in% curr.set, ]
    gen <- 1
    while (nrow(ped) > 0) {
        curr.set <- ped[!ped$sire %in% ped$id & !ped$dam %in% ped$id, ]
        for (i in 1:nrow(curr.set)) {
            tmp <- tbA[tbA$id1 %in% curr.set[i, 2:3] | tbA$id2 %in% curr.set[i, 2:3], ]
            tmp$a <- tmp$a/2
            tmp <- rbind(tmp, c(rep(curr.set[i, ]$id, 2), 1))
            tmp2 <- tmp[(tmp$id1 == curr.set[i, ]$sire & tmp$id2 == curr.set[i, ]$dam) | (tmp$id1 == curr.set[i, ]$dam & tmp$id2 == curr.set[i, ]$sire), ]
            if (nrow(tmp2) > 0) {
                tmp[nrow(tmp), ]$a <- tmp2$a + 1
                tmp = rbind(tmp, c(tmp2$id2, tmp2$id1, tmp2$a))
            }
            tmp[tmp$id1 %in% curr.set[i, 2:3], ]$id1 <- curr.set[i, ]$id
            if (nrow(tmp[tmp$id1 != curr.set[i, ]$id, ]) > 0)
                tmp[tmp$id1 != curr.set[i, ]$id, ]$id2 = curr.set[i, ]$id
            tmp[tmp$id1 > tmp$id2, 1:2] = tmp[tmp$id1 > tmp$id2, 2:1]
            tmp = aggregate(. ~ id1 + id2, data = tmp, sum)
            tbA = rbind(tbA, tmp)
        }
        ped = ped[!ped$id %in% curr.set$id, ]
        gen = gen + 1
    }
    tbA = tbA[order(tbA$id1, tbA$id2), ]

    sparseMatrix(i = tbA$id1, j = tbA$id2, x = tbA$a)

    message("Found ", gen, " generations")
    return(tbA)
}

