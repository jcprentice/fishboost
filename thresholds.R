library(data.table)
library(stringr)
library(ggplot2)

x <- fread("datasets/fb-final/data/scen-1-1-out/trace_combine.tsv")

# Extract just the genetic variances
x[, str_subset(names(x), "cov_G", negate = TRUE) := NULL]

# Get the factor associated with a top 10% / bottom 10%
factx <- function(x, p = 0.9) exp(2 * qnorm(p, mean = 0, sd = sqrt(x)))

di <- density(x$cov_G_ii)
ci <- rev(cumsum(rev(di$y))) / sum(di$y)
fi <- factx(di$x)

dt <- data.table(x = di$x, ci, fi)
mdt <- melt(dt, measure.vars = c("ci", "fi"))

ggplot(mdt) +
    geom_line(aes(x = x, y = value, colour = variable))
