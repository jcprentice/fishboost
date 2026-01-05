# read in Example 16 and investigate

{
    library(data.table)
    library(purrr)
    library(stringr)
    library(xml2)

    source("make_parameters.R")
    source("make_time_series.R")
    source("make_plots.R")
    source("get_ranks.R")
    source("get_roc.R")
}


dt <- read_xml("../sire21/examples/example16.xml") |>
    xml_child("datatable") |>
    xml_text() |>
    read.table(text = _) |>
    as.data.table()
setnames(dt, c("id", "group", "IC", "Tinf", "Tinc", "Tsym", "Tdeath", "inf", "sus", "tol"))

dt[, sdp := fifelse(is.na(group), "sire", "progeny")]
setcolorder(dt, c("id", "sdp"))
dt[, id := as.integer(str_remove(id, "Ind")) + 1]
dt[, IC := fifelse(IC == "S", 0L, 1L, NA_integer_)]
dt[, Tinf := as.numeric(fifelse(Tinf == "no", "0", Tinf, NA))]
dt[, Tinc := as.numeric(fifelse(Tinc == "no", "0", Tinc, NA))]
dt[, Tsym := as.numeric(fifelse(Tsym == "no", "0", Tsym, NA))]
dt[, Tdeath := as.numeric(fifelse(Tdeath == "no", "0", Tdeath, NA))]

params <- make_parameters(model_type = "SEIDR")

plt <- plot_model(dt, params)

estimated_BVs <- fread("../sire21/Output/estimated_bvs.csv")
posteriors    <- fread("../sire21/Output/posteriors.csv")
prec_accs     <- fread("../sire21/Output/pred_accs.csv")

true_vals <- c(beta = 0.05, tau_EI = 3, k_EI = 3, tau_EI = 4, k_EI = 3, tau_EI = 3, k_EI = 3,
               omega_gg = 1.0, omega_ff = 1.3, omega_rr = 0.5, omega_cor_gf = 0.4, omega_cor_gr = 0.1, omega_cor_fr = -0.2,
               sigma_gg = 0.7, sigma_ff = 1.0, sigma_rr = 0.9, sigma_cor_gf = 0.2, sigma_cor_gr = -0.1, sigma_cor_gf = -0.4)


# true_vals <- c(beta = 0.05, tau_EI = 3, k_EI = 3, tau_EI = 4, k_EI = 3, tau_EI = 3, k_EI = 3,
#                omega_gg = 1.0, omega_ff = 1.3,omega_cor_gf = 0.4,
#                sigma_gg = 0.7, sigma_ff = 1.0, sigma_cor_gf = 0.2)

posteriors[, true := true_vals]

posteriors[, .(parameter, true,
               mean = signif(mean, 3),
               rel_err = signif(abs(true-mean) / abs(true), 2))]

sires <- dt[sdp == "sire", id]
dt2 <- merge.data.table(dt[id %in% sires, .(id, sus, inf, tol)],
                        estimated_BVs[id %in% sires, .(id, est_sus = g_a, est_inf = f_a, est_tol = t_a)],
                        by = "id")

ranks <- dt2[, map(.SD, rank)]
# ranks[, `:=`(tol = seq_len(.N), est_tol = rank(sort(runif(.N)) + rnorm(.N, 0, 0.2)))]
plot_roc_curves(ranks, p = 0.8)

cor.test(ranks$sus, ranks$est_sus, method = "spearman")
cor.test(ranks$inf, ranks$est_inf, method = "spearman")
cor.test(ranks$tol, ranks$est_tol, method = "spearman")

