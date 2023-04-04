library(glue)
library(data.table)

x = readRDS("fb_data/fb_data12.rds")

x = x[sdp == "progeny"]
x[, c("sdp", "weight") := NULL]
x[trial == 1 & Trec == 104, Trec := NA]
x[trial == 2 & Trec == 160, Trec := NA]

val = "Tsym"
print(glue("Using {val} times:"))
cat("\n")

TP <- x[parasites == TRUE & !is.na(get(val)), .N]
FN <- x[parasites == TRUE & is.na(get(val)), .N]
FP <- x[parasites == FALSE & !is.na(get(val)), .N]
TN <- x[parasites == FALSE & is.na(get(val)), .N]

m1 <- matrix(c(TP, FN, FP, TN), 2, 2, byrow = TRUE,
             dimnames = list(parasites = c("yes", "no"),
                             died = c("yes", "no")))
print(m1)
cat("\n")
print(glue("Se = {se}, Sp = {sp}",
           se = round(100 * TP / (TP+FN), 1),
           sp = round(100 * TN / (TN+FP), 1)))
cat("\n")

# Repeat for donors only

donors = x[donor == 1]

TP <- donors[parasites == TRUE & !is.na(get(val)), .N]
FN <- donors[parasites == TRUE & is.na(get(val)), .N]
FP <- donors[parasites == FALSE & !is.na(get(val)), .N]
TN <- donors[parasites == FALSE & is.na(get(val)), .N]

m2 <- matrix(c(TP, FN, FP, TN), 2, 2, byrow = TRUE,
             dimnames = list(parasites = c("yes", "no"),
                             died = c("yes", "no")))
print(glue("Donors only:"))
cat("\n")
print(m2)
cat("\n")
print(glue("Se = {se}, Sp = {sp}",
           se = round(100 * TP / (TP+FN), 1),
           sp = round(100 * TN / (TN+FP), 1)))
