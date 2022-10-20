N <- 20
x <- data.table(id = (1:N) + 10)
x[, Tinf := runif(.N, 0, 10)]
x[, Tinc := Tinf + runif(.N, 0, 2)]
x[, Tsym := Tinc + runif(.N, 0, 2)]
x[, Trec := Tsym + runif(.N, 0, 5)]
x[]

epi_time <- 6
x[, status := factor("S", c("S", "E", "I", "D", "R"))]
x[Tinf < epi_time, status := "E"]
x[Tinc < epi_time, status := "I"]
x[Tsym < epi_time, status := "D"]
x[Trec < epi_time, status := "R"]
x[Tinf > epi_time, c("Tinf", "Tinc", "Tsym", "Trec") := 0]

print(table(x$status))

next_seidr_ni_event2 <- function(x, epi_time) {
    as.list(x[, .(.I,
                  a = fifelse(Tinc > epi_time, Tinc, Inf),
                  b = fifelse(Tsym > epi_time, Tsym, Inf),
                  c = fifelse(Trec > epi_time, Trec, Inf))
    ][, .(I, d = pmin(a, b, c))
    ][, .(t_next_event = min(d), id_next_event = which.min(d))])
}

next_seidr_ni_event(x, epi_time)
next_seidr_ni_event2(x, epi_time)
