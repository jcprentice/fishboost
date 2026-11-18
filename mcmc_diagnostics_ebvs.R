{
    library(data.table)
    library(purrr)
    library(stringr)
    library(coda)
    library(ggplot2)
}

mcmc_diagnostics_ebvs <- function(dataset = "fb-test",
                                  scen = 2,
                                  rep = 1,
                                  sdp = NULL) {
    if (FALSE) {
        dataset <- "fb-test"; scen <- 2; rep <- 1; sdp <- NULL
    }

    mod1 <- function(x, d) ((x - 1) %% d) + 1

    x <- readRDS(str_glue("datasets/{dataset}/data/scen-{scen}-{rep}-out/summary_inf.rds"))

    sdp1 <- if (!is.null(sdp)) sdp else c("sire", "dam", "progeny")
    y <- x$popn[sdp %in% sdp1, .SD, .SDcols = patterns("state|id|_g")]

    nchains <- 16L
    nstates <- last(x$popn$state) %/% nchains

    y[, chain := as.integer(ceiling(nchains * state / max(state)))]
    y[, state := mod1(state, nstates)]
    setcolorder(y, c("chain", "state"))

    z <- dcast(y, ... ~ id, value.var = str_subset(names(y), "_[ge]$"))
    zs <- split(z, by = "chain")

    nsample <- 2e6
    burnprop <- 0.2
    thinto <- nstates
    # w <- map(zs, ~ mcmc(.x[, .SD, .SDcols = patterns("_[ge]_")],))
    w <- map(zs, ~ mcmc(.x[, .SD, .SDcols = patterns("_[ge]_")],
                        start = nsample * burnprop,
                        thin = nsample / thinto))

    grd <- as.mcmc.list(w) |> gelman.diag(multivariate = FALSE) |> _$psrf[, 2]
    ess <- effectiveSize(w) |> round()

    list(grd = grd, ess = ess)
}

if (FALSE) {
    out_sires <- mcmc_diagnostics_ebvs("fb-test", 8, 1, "sire")
    out_all   <- mcmc_diagnostics_ebvs("fb-test", 8, 1)

    f <- function(x) {
        y <- c(max(x$grd), mean(x$grd), min(x$ess), mean(x$ess)) |> round(2)
        str_glue("GRD max  = {y[[1]]}",
                 "GRD mean = {y[[2]]}",
                 "ESS min  = {y[[3]]}",
                 "ESS mean = {y[[4]]}",
                 .sep = "\n", .trim = FALSE)
    }
    f(out_sires)
    f(out_all)
}
