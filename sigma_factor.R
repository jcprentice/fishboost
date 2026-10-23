# For a given variance s2, calculate how much more infectious
# the top 5% (say) are than the lowest 5%

# Algebraic result
sigma2f <- function(s2 = 1, p = 0.05) {
    s <- sqrt(s2)
    z <- qnorm(1 - q)
    t <- exp(- z^2 / 2) / (sqrt(2 * pi) * p)
    exp(2 * s * t)
}

# Numerical result
check_s2f <- function(s2 = 1, p = 0.05) {
    s <- sqrt(s2)
    # Sample some values with sigma
    x <- rnorm(1e6, 0, s)
    # Note that the distribution is symmetrical,
    # so only require top q, as top_q = - bottom_q
    top_q <- quantile(x, 1 - p)
    m <- mean(x[x > top_q])
    exp(2 * m)
}

s2 <- 0.5
q <- 0.05

# These should be roughly similar
sigma2f(s2, q)
check_s2f(s2, q)


top20 <- function(s2, N = 1e3) {
    s <- sqrt(s2)

    # Threshold for top 20%
    q80 <- qlnorm(0.8, 0, s)

    # Equally spaced points excluding (0, Inf)
    qs <- seq(0, 1, 1/N)[-c(1, N+1)]
    x <- qlnorm(qs, 0, s)

    # Proportion above q80
    sum(x[x > q80]) / sum(x)
}


top20(2.83497, 1e6)
