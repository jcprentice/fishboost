get_R0 <- function(pop) {
    x <- table(pop$generation)

    R0 <- if (length(x) > 1) x[[2]] / x[[1]] else 0

    message("R0 estimate: ", signif(R0, 3))

    R0
}

