x = readRDS("fb_data/fb_12.rds")

x = x[sdp == "progeny"]

trial_size = x[, .N, trial][, N]
nprogeny = sum(trial_size)
ndonors = x[donor == 1, .N]

x[, parasites := factor(fifelse(parasites, "parasites", "no parasites"), c("parasites", "no parasites"))]
tmax = c(104, 160)
x[, died := factor(fifelse(Tdeath %notin% tmax[trial], "died", "survived"), c("died", "survived"))]
x[, symptoms := factor(fifelse(!is.na(Tsym), "symptoms", "no symptoms"), c("symptoms", "no symptoms"))]

message("parasites / died")
print(x[, .(.N, pc=round(.N * 100 / nprogeny, 1)), .(parasites, died)])

message("parasites / symptoms")
print(x[, .(.N, pc = round(.N * 100 / nprogeny, 1)), .(parasites, symptoms)])

message("parasites / died, donors only")
print(x[donor == 1, .(.N, pc = round(.N * 100 / ndonors, 1)), .(parasites, died)])

message("parasites / symptoms, donors only")
print(x[donor == 1, .(.N, pc = round(.N * 100 / ndonors, 1)), .(parasites, symptoms)])


message("parasites / died / trial")
print(x[, .(.N, pc = round(.N * 100 / trial_size[trial], 1)), .(parasites, died, trial)])

message("parasites / symptoms / trial")
print(x[, .(.N, pc = round(.N * 100 / trial_size[trial], 1)), .(parasites, symptoms, trial)])

