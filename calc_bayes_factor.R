library(data.table)
library(stringr)
library(purrr)

source("utils.R")

calc_BF <- function(data, param = "cov_G_ii", params, pt = 0, pc = 1) {
    # pt: the fixed point to compare models (usually 0, but could take other
    # numbers, e.g. 1 for the shape of a gamma distribution)
    #
    # pc: the percentage of the distribution around the fixed point
    
    pc <- clamp(pc / 100, 0, 1)
    prior_rng <- params$priors[parameter == param, c(val1, val2)]
    width <- abs(diff(prior_rng))
    
    # Can't calculate density for a fixed point (var = 0), so just return NA
    if (var(data[[param]]) == 0) {
        return(NA_real_)
    }
    
    dens <- density(data[[param]],
                    from = min(pt, prior_rng[[1]]),
                    to   = max(pt, prior_rng[[2]]),
                    bw = "SJ", adjust = 1, cut = 0)
    dd <- with(dens, data.table(x, y))
    area <- approx(dd$x, dd$y, pt)$y * width
    signif(1 / area, 3)
}


calc_bayes_factor <- function(dataset = "fb-final", scen = 1, pc = 1) {
    # dataset <- "fb-final"; scen <- 1; pc <- 1
    
    res <- readRDS(str_glue("datasets/{dataset}/results/scen-{scen}-1.rds"))
    params <- res$params
    data <- fread(str_glue("datasets/{dataset}/data/scen-{scen}-1-out/trace_combine.tsv"))
    
    data[, str_subset(names(.SD), "Group|state") := NULL]
    param_names <- names(data)
    
    map_dbl(param_names, \(param) {
        # param <- "cov_G_ii"
        pt <- if (str_detect(param, "shape")) 1 else 0
        BF <- calc_BF(data, param, params, pt, pc)
                       
        # every sqrt(10) is worth 1 star
        stars <- str_dup("*", min(abs(log10(BF)), 2) %/% 0.5)
        
        message(str_glue("{param} = {BF} {stars}"))
        BF
    }) |> set_names(param_names)
}

# dataset <- "fb-weight-grm"; scen <- 2
dataset <- "fb-final"; scen <- 1

BFs <- calc_bayes_factor(dataset, scen)

