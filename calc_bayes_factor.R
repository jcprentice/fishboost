{
    library(data.table)
    library(stringr)
    library(purrr)

    safe_source("utils.R")
}

calc_BF <- function(par = "cov_G_ii", data, params, pt = 0, qnt = 0.01) {
    # pt: the fixed point to compare models (usually 0, but could take other
    # numbers, e.g. 1 for the shape of a gamma distribution)
    #
    # qnt: the quantile of the distribution around the fixed point

    # Can't calculate density for a fixed point (var = 0), so just return NA
    if (par %notin% params$priors$parameter || var(data[[par]]) == 0) {
        return(NA_real_)
    }

    # Get prior parameters
    type <- params$priors[parameter == par, type]
    A <- params$priors[parameter == par, val1]
    B <- params$priors[parameter == par, val2]

    pt <- clamp(pt, A, B)
    c(A = A, pt = pt, B = B)

    prior_dens <- if (type == "default") {
        if (str_starts(par, "cov_")) {
            # Half normal(0, 2)
            dnorm(pt, 0, 2) * 2
        } else {
            # LKJ(1.2)
            dbeta(pt / 2 + 0.5, 1.2, 1.2)
        }
    } else if (type == "inverse") {
        1 / (pt * log(B / A))
    } else if (type == "uniform") {
        dunif(pt, A, B)
    }
    # prior_dens

    # post_dens <- mean(between(data[[par]], a, b))


    dens <- density(data[[par]],
                    adjust = 1, # cut = 0, # bw = "SJ",
                    from = A, to = B)
    dd <- with(dens, data.table(x, y))
    post_dens <- approx(dd$x, dd$y, pt)$y

    BF <- post_dens / prior_dens
    BF
}


calc_bayes_factor <- function(dataset = "fb-final", scen = 1, qnt = 0.01) {
    if (FALSE) {
        dataset <- "fb-test"; scen <- 2; qnt <- 0.01
    }

    res <- readRDS(str_glue("datasets/{dataset}/results/scen-{scen}-1.rds"))
    params <- res$params
    data <- fread(str_glue("datasets/{dataset}/data/scen-{scen}-1-out/trace_combine.tsv"))

    data[, str_subset(names(.SD), "Group|state") := NULL]
    pars <- names(data) |>
        str_subset("^G_", negate = TRUE)

    BFs <- map(pars, \(par) {
        if (FALSE) {
            i <- 1; par <- pars[[i]]
        }

        pt <- if (str_detect(par, "[LDR]P|shape|infrat")) {
            1
        } else if (str_detect(par, "beta")) {
            0.5
        } else {
            0
        }

        BF <- calc_BF(par, data, params, pt, qnt) |> signif(3)

        # every sqrt(10) is worth 1 star
        stars <- if (is.na(BF)) {
            NA_character_
        } else if (BF > 1) {
            str_dup("!", min(log10(BF), 1.5) %/% 0.5)
        } else {
            str_dup("*", min(-log10(BF), 1.5) %/% 0.5)
        }

        # message(str_glue("{par} != {pt}:\t{BF}\t{stars}"))
        list(par = par, pt = pt, BF = BF, stars = stars)
    }) |>
        rbindlist()

    BFs[!is.na(BF), .(par, pt, BF = signif(BF, 3), logBF = round(log10(BF), 2), stars)]
}

if (FALSE) {
    BFs_fb_test_s1 <- calc_bayes_factor("fb-test", 1)
    BFs_fb_test_s2 <- calc_bayes_factor("fb-test", 2)
    BFs_fb_test_s8 <- calc_bayes_factor("fb-test", 8)
}

