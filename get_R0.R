get_R0 <- function(pop) {
    x <- table(pop$generation)
    if (length(x) > 1)
        R0 <- x[[2]] / x[[1]]
    else
        R0 <- 0
    message("R0 estimate: ", signif(R0, 3))

    R0
}

