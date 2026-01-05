# Algebraic result
sigma2f <- function(sigma = 1, q = 0.95) {
    z <- qnorm(q)
    t <- exp(-z^2/2) / (sqrt(2*pi) * (1-q))
    exp(2 * sigma * t)
}

# Numerical result
check_s2f <- function(sigma = 1, q = 0.95) {
   # Sample some values with sigma 
   x <- rnorm(1e6, 0, sd = sigma)
   # Note that the distribution is symmetrical,
   # so only require top q, as top_q = - bottom_q
   top_q <- quantile(x, q)
   m <- mean(x[x > top_q])
   # exp(top_m) / exp(bottom_m)
   # = exp(top_m - bottom_m)
   # = exp(2 * top_m)
   exp(2 * m)
}

sigma <- 0.5
q <- 0.95

sigma2f(sigma, q)
check_s2f(sigma, q)

