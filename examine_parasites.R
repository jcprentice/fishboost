x = readRDS("fb_data/fb_12.rds")

x = x[sdp == "progeny"]

trial_size = x[, .N, trial][, N]
nprogeny = sum(trial_size)
ndonors = x[donor == 1, .N]

x[, parasites := factor(fifelse(parasites, "parasites", "no parasites"), c("parasites", "no parasites"))]
tmax = c(104, 160)
x[, died := factor(fifelse(Tdeath %notin% tmax[trial], "died", "survived"), c("died", "survived"))]
x[, signs := factor(fifelse(!is.na(Tsign), "signs", "no signs"), c("signs", "no signs"))]

message("parasites / died")
print(x[, .(.N, pc = round(.N * 100 / nprogeny, 1)), .(parasites, died)])

message("parasites / signs")
print(x[, .(.N, pc = round(.N * 100 / nprogeny, 1)), .(parasites, signs)])

message("parasites / died, donors only")
print(x[donor == 1, .(.N, pc = round(.N * 100 / ndonors, 1)), .(parasites, died)])

message("parasites / signs, donors only")
print(x[donor == 1, .(.N, pc = round(.N * 100 / ndonors, 1)), .(parasites, signs)])


message("parasites / died / trial")
print(x[, .(.N, pc = round(.N * 100 / trial_size[trial], 1)), .(parasites, died, trial)])

message("parasites / signs / trial")
print(x[, .(.N, pc = round(.N * 100 / trial_size[trial], 1)), .(parasites, signs, trial)])

