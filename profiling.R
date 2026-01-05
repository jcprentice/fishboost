tmp <- tempfile()
Rprof(tmp, interval = 0.1)
popn <- simulate_epidemic(popn, params)
Rprof(NULL)

