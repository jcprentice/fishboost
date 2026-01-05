library(data.table)
library(stringr)

widen_priors <- function(priors) {
    priors[str_starts(parameter, "cov_"),
           `:=`(val1 = 0,
                val2 = pmax(val2, ceiling(true_val + 1), na.rm = TRUE))]
    priors[str_starts(parameter, "r_"),
           `:=`(val1 = -0.95, val2 = 0.95)]
    priors[str_detect(parameter, "period|beta|sigma"),
           `:=`(val1 = 0,
                val2 = pmax(val2, true_val * 1.2, na.rm = TRUE))]
    priors[str_ends(parameter, "shape"),
           `:=`(val1 = 0.5, val2 = 5)]
    priors[str_starts(parameter, "trial|donor|weight"),
           `:=`(val1 = pmin(val1, floor(true_val - 2), na.rm = TRUE),
                val2 = pmax(val2, ceiling(true_val + 2), na.rm = TRUE))]
    
    priors
}
