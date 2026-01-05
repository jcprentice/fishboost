library(data.table)
library(fitdistrplus)
library(purrr)

get_RPs <- function() {
    km_data <- readRDS("fb-final-old/results/km.rds")
    x <- km_data$km_data[[1]]$data[src == "sim", .(trial, donor, Tsym, Tdeath)]
    
    # Storage
    dt <- list(Trial = 1:2,
               Period = "RP",
               Subset = c("All", "Seeders", "Contacts"),
               mean = 0, rate = 0, shape = 0) |>
        rev() |> expand.grid(stringsAsFactors = FALSE) |> rev() |> setDT()
    vals <-  c("mean", "rate", "shape")
    
    
    # Cut out values where both Tsym and Tdeath2 are NA (need at least 1 value)
    x <- x[!is.na(Tsym) | !is.na(Tdeath)]
    
    # Trial 1
    f1 <- x[trial == 1 & !is.na(Tsym) &
                 !is.na(Tdeath) & Tdeath > Tsym, Tdeath - Tsym] |>
        fitdist("gamma")
    f1$estimate[["mean"]] = f1$estimate[["shape"]] / f1$estimate[["rate"]]
    message("Trial 1:")
    print(round(f1$estimate[vals], 2))
    set(dt, 1L, vals, as.list(f1$estimate[vals]))
    
    # Trial 1 seeders
    f1s <- x[trial == 1 & donor == 1 & !is.na(Tsym) &
                 !is.na(Tdeath) & Tdeath > Tsym, Tdeath - Tsym] |>
        fitdist("gamma")
    f1s$estimate[["mean"]] = f1s$estimate[["shape"]] / f1s$estimate[["rate"]]
    message("Trial 1 seeders:")
    print(round(f1s$estimate[vals], 2))
    set(dt, 2L, vals, as.list(f1s$estimate[vals]))
    
    # Trial 1 contact
    f1c <- x[trial == 1 & donor == 0 & !is.na(Tsym) &
                 !is.na(Tdeath) & Tdeath > Tsym, Tdeath - Tsym] |>
        fitdist("gamma")
    f1c$estimate[["mean"]] = f1c$estimate[["shape"]] / f1c$estimate[["rate"]]
    message("Trial 1 contact:")
    print(round(f1c$estimate[vals], 2))
    set(dt, 3L, vals, as.list(f1c$estimate[vals]))
    
    # Trial 2
    f2 <- x[trial == 2 & !is.na(Tsym) &
                 !is.na(Tdeath) & Tdeath > Tsym, Tdeath - Tsym] |>
        fitdist("gamma")
    f2$estimate[["mean"]] = f2$estimate[["shape"]] / f2$estimate[["rate"]]
    message("Trial 2:")
    print(round(f2$estimate[vals], 2))
    set(dt, 4L, vals, as.list(f2$estimate[vals]))
    
    # Trial 2 seeders
    f2s <- x[trial == 2 & donor == 1 & !is.na(Tsym) &
                 !is.na(Tdeath) & Tdeath > Tsym, Tdeath - Tsym] |>
        fitdist("gamma")
    f2s$estimate[["mean"]] = f2s$estimate[["shape"]] / f2s$estimate[["rate"]]
    message("Trial 2 seeders:")
    print(round(f2s$estimate[vals], 2))
    set(dt, 5L, vals, as.list(f2s$estimate[vals]))
    
    # Trial 2 contact
    f2c <- x[trial == 2 & donor == 0 & !is.na(Tsym) &
                 !is.na(Tdeath) & Tdeath > Tsym, Tdeath - Tsym] |>
        fitdist("gamma")
    f2c$estimate[["mean"]] = f2c$estimate[["shape"]] / f2c$estimate[["rate"]]
    message("Trial 2 contact:")
    print(round(f2c$estimate[vals], 2))
    set(dt, 6L, vals, as.list(f2c$estimate[vals]))
    
    mget(c("dt", "f1", "f2", "f1s", "f1c", "f2s", "f2c"))
}

get_RPs_fb <- function(fix_seeders = c(10, 80)) {
    # Process data
    fb_data <- readRDS("fb_data/fb_12.rds")
    x <- fb_data[sdp == "progeny", .(trial, donor, Tsym, Tdeath, parasites)]
    
    x[Tdeath == c(104, 160)[trial], Tdeath := NA]
    if (any(fix_seeders > 0)) {
        message("Reclassifying seeders")
        x[donor == 1 & Tsym > fix_seeders[trial], donor := 0]
    }
    
    # Storage
    dt <- list(Trial = 1:2,
               Period = "RP",
               Subset = c("All", "Seeders", "Contacts"),
               mean = 0, rate = 0, shape = 0) |>
        rev() |> expand.grid(stringsAsFactors = FALSE) |> rev() |> setDT()
    
    vals <-  c("mean", "rate", "shape")
    
    # Cut out values where both Tsym and Tdeath2 are NA (need at least 1 value)
    x <- x[!is.na(Tsym) | !is.na(Tdeath)]
    
    x[, `:=`(left = fcase(!is.na(Tsym) & !is.na(Tdeath), pmax(Tdeath - Tsym - 1, 0),
                          !is.na(Tsym) & is.na(Tdeath), c(104, 160)[trial] - Tsym,
                          is.na(Tsym) & !is.na(Tdeath), 0,
                          default = NA),
             right = fcase(!is.na(Tsym) & !is.na(Tdeath), pmax(Tdeath - Tsym + 1, 0),
                           !is.na(Tsym) & is.na(Tdeath), NA,
                           is.na(Tsym) & !is.na(Tdeath), Tdeath,
                           default = NA))]
    
    # Trial 1
    tmp <- x[trial == 1 & !is.na(Tsym) &
                 !is.na(Tdeath) & Tdeath > Tsym, Tdeath - Tsym] |>
        fitdist("gamma")
    f1 <- fitdistcens(x[trial == 1, .(left, right)], "gamma",
                       start = list(shape = tmp$estimate[["shape"]],
                                    rate = tmp$estimate[["rate"]]))
    f1$estimate[["mean"]] = f1$estimate[["shape"]] / f1$estimate[["rate"]]
    message("Trial 1:")
    print(round(f1$estimate[vals], 2))
    set(dt, 1L, vals, as.list(f1$estimate[vals]))
    
    # Trial 1 seeders
    tmp <- x[trial == 1 & donor == 1 & !is.na(Tsym) &
                 !is.na(Tdeath) & Tdeath > Tsym, Tdeath - Tsym] |>
        fitdist("gamma")
    f1s <- fitdistcens(x[trial == 1 & donor == 1, .(left, right)], "gamma",
                       start = list(shape = tmp$estimate[["shape"]],
                                    rate = tmp$estimate[["rate"]]))
    f1s$estimate[["mean"]] = f1s$estimate[["shape"]] / f1s$estimate[["rate"]]
    message("Trial 1 seeders:")
    print(round(f1s$estimate[vals], 2))
    set(dt, 2L, vals, as.list(f1s$estimate[vals]))
    
    # Trial 1 contact
    tmp <- x[trial == 1 & donor == 0 & !is.na(Tsym) &
                 !is.na(Tdeath) & Tdeath > Tsym, Tdeath - Tsym] |>
        fitdist("gamma")
    f1c <- fitdistcens(x[trial == 1 & donor == 0, .(left, right)], "gamma",
                       start = list(shape = tmp$estimate[["shape"]],
                                    rate = tmp$estimate[["rate"]]))
    f1c$estimate[["mean"]] = f1c$estimate[["shape"]] / f1c$estimate[["rate"]]
    message("Trial 1 contact:")
    print(round(f1c$estimate[vals], 2))
    set(dt, 3L, vals, as.list(f1c$estimate[vals]))
    
    # Trial 2
    tmp <- x[trial == 2 & !is.na(Tsym) &
                 !is.na(Tdeath) & Tdeath > Tsym, Tdeath - Tsym] |>
        fitdist("gamma")
    f2 <- fitdistcens(x[trial == 2, .(left, right)], "gamma",
                       start = list(shape = tmp$estimate[["shape"]],
                                    rate = tmp$estimate[["rate"]]))
    f2$estimate[["mean"]] = f2$estimate[["shape"]] / f2$estimate[["rate"]]
    message("Trial 2:")
    print(round(f2$estimate[vals], 2))
    set(dt, 4L, vals, as.list(f2$estimate[vals]))
    
    # Trial 2 seeders
    tmp <- x[trial == 2 & donor == 1 & !is.na(Tsym) &
                 !is.na(Tdeath) & Tdeath > Tsym, Tdeath - Tsym] |>
        fitdist("gamma")
    f2s <- fitdistcens(x[trial == 2 & donor == 1, .(left, right)], "gamma",
                       start = list(shape = tmp$estimate[["shape"]],
                                    rate = tmp$estimate[["rate"]]))
    f2s$estimate[["mean"]] = f2s$estimate[["shape"]] / f2s$estimate[["rate"]]
    message("Trial 2 seeders:")
    print(round(f2s$estimate[vals], 2))
    set(dt, 5L, vals, as.list(f2s$estimate[vals]))
    
    # Trial 2 contact
    tmp <- x[trial == 2 & donor == 0 & !is.na(Tsym) &
                 !is.na(Tdeath) & Tdeath > Tsym, Tdeath - Tsym] |>
        fitdist("gamma")
    f2c <- fitdistcens(x[trial == 2 & donor == 0, .(left, right)], "gamma",
                       start = list(shape = tmp$estimate[["shape"]],
                                    rate = tmp$estimate[["rate"]]))
    f2c$estimate[["mean"]] = f2c$estimate[["shape"]] / f2c$estimate[["rate"]]
    message("Trial 2 contact:")
    print(round(f2c$estimate[vals], 2))
    set(dt, 6L, vals, as.list(f2c$estimate[vals]))
    
    mget(c("dt", "f1", "f2", "f1s", "f1c", "f2s", "f2c"))
}


get_DPs <- function() {
    km_data <- readRDS("fb-final-old/results/km.rds")
    x <- km_data$km_data[[1]]$data[src == "sim", .(trial, donor, Tsym)]
    
    # Storage
    dt <- list(Trial = 1:2,
               Period = "Tsym",
               Subset = c("All", "Seeders", "Contacts"),
               mean = 0, rate = 0, shape = 0) |>
        rev() |> expand.grid(stringsAsFactors = FALSE) |> rev() |> setDT()
    
    vals <-  c("mean", "rate", "shape")
    
    # Cut out values where Tsym is NA
    x <- x[!is.na(Tsym)]
    
    # Trial 1
    f1 <- x[trial == 1 & !is.na(Tsym), Tsym] |>
        fitdist("gamma")
    f1$estimate[["mean"]] = f1$estimate[["shape"]] / f1$estimate[["rate"]]
    message("Trial 2 contact:")
    print(round(f1$estimate[vals], 2))
    set(dt, 1L, vals, as.list(f1$estimate[vals]))
    
    # Trial 1 seeders
    f1s <- x[trial == 1 & donor == 1 & !is.na(Tsym), Tsym] |>
        fitdist("gamma")
    f1s$estimate[["mean"]] = f1s$estimate[["shape"]] / f1s$estimate[["rate"]]
    message("Trial 1 seeders:")
    print(round(f1s$estimate[vals], 2))
    set(dt, 2L, vals, as.list(f1s$estimate[vals]))
    
    # Trial 1 contact
    f1c <- x[trial == 1 & donor == 0 & !is.na(Tsym), Tsym] |>
        fitdist("gamma")
    f1c$estimate[["mean"]] = f1c$estimate[["shape"]] / f1c$estimate[["rate"]]
    message("Trial 1 contact:")
    print(round(f1c$estimate[vals], 2))
    set(dt, 3L, vals, as.list(f1c$estimate[vals]))
    
    # Trial 2
    f2 <- x[trial == 2 & !is.na(Tsym), Tsym] |>
        fitdist("gamma")
    f2$estimate[["mean"]] = f2$estimate[["shape"]] / f2$estimate[["rate"]]
    message("Trial 2 contact:")
    print(round(f2$estimate[vals], 2))
    set(dt, 4L, vals, as.list(f2$estimate[vals]))
    
    # Trial 2 seeders
    f2s <- x[trial == 2 & donor == 1 & !is.na(Tsym), Tsym] |>
        fitdist("gamma")
    f2s$estimate[["mean"]] = f2s$estimate[["shape"]] / f2s$estimate[["rate"]]
    message("Trial 2 seeders:")
    print(round(f2s$estimate[vals], 2))
    set(dt, 5L, vals, as.list(f2s$estimate[vals]))
    
    # Trial 2 contact
    f2c <- x[trial == 2 & donor == 0 & !is.na(Tsym), Tsym] |>
        fitdist("gamma")
    f2c$estimate[["mean"]] = f2c$estimate[["shape"]] / f2c$estimate[["rate"]]
    message("Trial 2 contact:")
    print(round(f2c$estimate[vals], 2))
    set(dt, 6L, vals, as.list(f2c$estimate[vals]))
    
    mget(c("dt", "f1", "f2", "f1s", "f1c", "f2s", "f2c"))
}

get_DPs_fb <- function(fix_seeders = c(10, 80)) {
    # Process data
    fb_data <- readRDS("fb_data/fb_12.rds")
    x <- fb_data[sdp == "progeny", .(trial, donor, Tsym, Tdeath)]
    
    x[Tdeath == c(104, 160)[trial], Tdeath := NA]
    if (any(fix_seeders > 0)) {
        message("Reclassifying seeders")
        x[donor == 1 & Tsym > fix_seeders[trial], donor := 0]
    }
    
    # Storage
    dt <- list(Trial = 1:2,
               Period = "Tsym",
               Subset = c("All", "Seeders", "Contacts"),
               mean = 0, rate = 0, shape = 0) |>
        rev() |> expand.grid(stringsAsFactors = FALSE) |> rev() |> setDT()
    
    vals <-  c("mean", "rate", "shape")
    
    # Cut out values where both Tsym and Tdeath2 are NA (need at least 1 value)
    x <- x[!is.na(Tsym) | !is.na(Tdeath)]
    
    x[, `:=`(left = fcase(!is.na(Tsym), pmax(Tsym - 1, 0),
                          is.na(Tsym) & !is.na(Tdeath), 0,
                          default = NA),
             right = fcase(!is.na(Tsym), Tsym,
                           is.na(Tsym) & !is.na(Tdeath), Tdeath,
                           default = NA))]
    
    # Trial 1
    tmp <- x[trial == 1 & !is.na(Tsym), Tsym] |>
        fitdist("gamma")
    f1 <- fitdistcens(x[trial == 1, .(left, right)], "gamma",
                       start = list(shape = tmp$estimate[["shape"]],
                                    rate = tmp$estimate[["rate"]]))
    f1$estimate[["mean"]] = f1$estimate[["shape"]] / f1$estimate[["rate"]]
    message("Trial 2 contact:")
    print(round(f1$estimate[vals], 2))
    set(dt, 1L, vals, as.list(f1$estimate[vals]))
    
    # Trial 1 seeders
    tmp <- x[trial == 1 & donor == 1 & !is.na(Tsym), Tsym] |>
        fitdist("gamma")
    f1s <- fitdistcens(x[trial == 1 & donor == 1, .(left, right)], "gamma",
                       start = list(shape = tmp$estimate[["shape"]],
                                    rate = tmp$estimate[["rate"]]))
    f1s$estimate[["mean"]] = f1s$estimate[["shape"]] / f1s$estimate[["rate"]]
    message("Trial 1 seeders:")
    print(round(f1s$estimate[vals], 2))
    set(dt, 2L, vals, as.list(f1s$estimate[vals]))
    
    # Trial 1 contact
    tmp <- x[trial == 1 & donor == 0 & !is.na(Tsym), Tsym] |>
        fitdist("gamma")
    f1c <- fitdistcens(x[trial == 1 & donor == 0, .(left, right)], "gamma",
                       start = list(shape = tmp$estimate[["shape"]],
                                    rate = tmp$estimate[["rate"]]))
    f1c$estimate[["mean"]] = f1c$estimate[["shape"]] / f1c$estimate[["rate"]]
    message("Trial 1 contact:")
    print(round(f1c$estimate[vals], 2))
    set(dt, 3L, vals, as.list(f1c$estimate[vals]))
    
    # Trial 2
    tmp <- x[trial == 2 & !is.na(Tsym), Tsym] |>
        fitdist("gamma")
    f2 <- fitdistcens(x[trial == 2, .(left, right)], "gamma",
                       start = list(shape = tmp$estimate[["shape"]],
                                    rate = tmp$estimate[["rate"]]))
    f2$estimate[["mean"]] = f2$estimate[["shape"]] / f2$estimate[["rate"]]
    message("Trial 2 contact:")
    print(round(f2$estimate[vals], 2))
    set(dt, 4L, vals, as.list(f2$estimate[vals]))
    
    # Trial 2 seeders
    tmp <- x[trial == 2 & donor == 1 & !is.na(Tsym), Tsym] |>
        fitdist("gamma")
    f2s <- fitdistcens(x[trial == 2 & donor == 1, .(left, right)], "gamma",
                       start = list(shape = tmp$estimate[["shape"]],
                                    rate = tmp$estimate[["rate"]]))
    f2s$estimate[["mean"]] = f2s$estimate[["shape"]] / f2s$estimate[["rate"]]
    message("Trial 2 seeders:")
    print(round(f2s$estimate[vals], 2))
    set(dt, 5L, vals, as.list(f2s$estimate[vals]))
    
    # Trial 2 contact
    tmp <- x[trial == 2 & donor == 0 & !is.na(Tsym), Tsym] |>
        fitdist("gamma")
    f2c <- fitdistcens(x[trial == 2 & donor == 0, .(left, right)], "gamma",
                       start = list(shape = tmp$estimate[["shape"]],
                                    rate = tmp$estimate[["rate"]]))
    f2c$estimate[["mean"]] = f2c$estimate[["shape"]] / f2c$estimate[["rate"]]
    message("Trial 2 contact:")
    print(round(f2c$estimate[vals], 2))
    set(dt, 6L, vals, as.list(f2c$estimate[vals]))
    
    mget(c("dt", "f1", "f2", "f1s", "f1c", "f2s", "f2c"))
}

DPs <- get_DPs()
RPs <- get_RPs()

DPsfb1 <- get_DPs_fb(fix_seeders = c(0, 0))
DPsfb2 <- get_DPs_fb(fix_seeders = c(10, 80))

RPsfb1 <- get_RPs_fb(fix_seeders = c(0, 0))
RPsfb2 <- get_RPs_fb(fix_seeders = c(10, 80))

Ps <- rbind(DPs$dt, RPs$dt)
Psfb1 <- rbind(DPsfb1$dt, RPsfb1$dt)
Psfb2 <- rbind(DPsfb2$dt, RPsfb2$dt)

setorder(Ps, "Trial")
setorder(Psfb1, "Trial")
setorder(Psfb2, "Trial")

cols <- c("mean", "rate", "shape")
Ps[, (cols) := map(.SD, round, 2), .SDcols = cols]
Psfb1[, (cols) := map(.SD, round, 2), .SDcols = cols]
Psfb2[, (cols) := map(.SD, round, 2), .SDcols = cols]
