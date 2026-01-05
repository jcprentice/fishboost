get_R0 <- function(popn) {
    x <- table(popn$generation)

    R0 <- if (length(x) > 1) x[[2]] / x[[1]] else 0

    message("R0 estimate: ", signif(R0, 3))

    R0
}

