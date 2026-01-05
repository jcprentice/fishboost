library(ggplot2)
library(stringr)
library(data.table)

# Given a multiplier, find the variance s^2
binsearch <- function(f, q = 0.9, eps = 1e-6) {
    L <- 0
    H <- 10
    
    # Limit iterations to this value
    itns <- 30
    while (itns > 0) {
        itns <- itns - 1
        M <- (L + H) / 2
        prop <- exp(2 * qnorm(q, 0, sqrt(M)))
        
        if (abs(prop - f) < eps * abs(f)) return(M)
        else if (prop < f) L <- M
        else if (prop > f) H <- M
    }
    NA_real_
}

# sigma = 0.5, mean of top 10% = 0.878, factor = 5.79x

topq <- function(q = 0.9, mu = 0, sigma = 1) {
    z <- qnorm(q)
    t <- exp(-z^2 / 2) / (sqrt(2 * pi) * (1 - q))
    sigma * t + mu
}



sigma2f <- function(sigma = 0.5, mu = 0, q = 0.9) {
    exp(2 * topq(q, mu, sigma))
}


f2sigma <- function(f, q = 0.9) {
    z <- qnorm(q)
    log(f) * (1 - q) * sqrt(2 * pi) / 2 * exp(z^2 / 2)
}


get_mult <- function(dataset = "fb-final", scen = 1) {
    
    # Load data and trim unwanted columns
    x <- fread(str_glue("datasets/{dataset}/data/scen-{scen}-1-out/trace_combine.tsv"))
    x[, str_subset(names(x), "Group|state") := NULL]
    x[, R0 := cov_G_ss + cov_G_ii + cov_G_tt]
    
    # Ideally I want to solve the following equation for sd given f, but I
    # can't, so I do a search-based solve instead
    
    # f = exp(2 * qnorm(0.95, 0, sd))
    # log(f) / 2 = qnorm(0.95, 0, sd)
    # sd = ?
    
    
    # Build table of multipliers and the probability
    # dt2 <- data.table(multiplier = c(1.2, 1.5, 2, 5, 10), p = 0)
    dt2 <- data.table(multiplier = seq(1.1, 10, by = 0.1),
                      p1 = 0, p2 = 0)
    q1 <- 0.95
    q2 <- 0.95
    walk(seq_len(nrow(dt2)), \(i) {
        m <- dt2$multiplier[i]
        s1 <- binsearch(m, q1)
        p1 <- mean(x$R0 > s1)
        s2 <- f2sigma(m, q2)
        p2 <- mean(x$R0 > s2)
        set(dt2, i, c("p1", "p2"), list(p1, p2))
    })
    
    dt2
}

# dataset <- "fb-parasites4"; scen <- 3
dataset <- "fb-fes3"; scen <- 1

m <- get_mult(dataset, scen)
m[c(2, 5, 10, 40, 90), .(multiplier, p1 = round(p1, 3), p2 = round(p2, 3))]

m2 <- melt(m, id.vars = "multiplier")


plt <- ggplot(m2) +
    geom_line(aes(x = multiplier, y = value, colour = variable),
              linewidth = 1) +
    labs(x = "Factor", y = "p >= Factor") +
    lims(y = c(0, 1))
ggsave("multiplier-R0.png", plt, width = 8, height = 4)

