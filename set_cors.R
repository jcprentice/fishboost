library(stringr)

set_cors <- function(X, cors) {
    sit <- c("susceptibility", "infectivity", "tolerance")
    nx <- rownames(X)
    
    for (i in seq_len(nrow(X))) for (j in seq(i, ncol(X))) {
        if (i == j) {
            next()
        } else if (!nx[i] %in% sit || !nx[j] %in% sit) {
            X[j, i] <- X[i, j] <- 0
        } else {
            X[j, i] <- X[i, j] <- cors
        }
    }
    X
}
