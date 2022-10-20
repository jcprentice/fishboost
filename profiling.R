tmp <- tempfile()
Rprof(tmp, interval = 0.1)
system.time(pop <- simulate_epidemic(traits, params))
Rprof(NULL)
