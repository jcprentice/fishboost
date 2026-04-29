x <- readRDS("datasets/fb-test/data/scen-7-1-out/summary_inf.rds")$popn

sires_s <- x[sdp == "sire", .(id, trait = sus_g)]
sires_i <- x[sdp == "sire", .(id, trait = inf_g)]
sires_t <- x[sdp == "sire", .(id, trait = tol_g)]

ids_s <- sires_s[, median(trait), id][order(V1), id[c(1, .N)]]
ids_i <- sires_i[, median(trait), id][order(V1), id[c(1, .N)]]
ids_t <- sires_t[, median(trait), id][order(V1), id[c(1, .N)]]

ratios_s <- sires_s[id == ids_s[2], exp(trait)] / sires_s[id == ids_s[1], exp(trait)]
ratios_i <- sires_i[id == ids_i[2], exp(trait)] / sires_i[id == ids_i[1], exp(trait)]
ratios_t <- sires_t[id == ids_t[2], exp(trait)] / sires_t[id == ids_t[1], exp(trait)]

N <- length(ratios_s)

traits <- c("sus", "inf", "tol")
dt <- data.table(trait = rep(traits, each = N),
                 ratio = c(ratios_s, ratios_i, ratios_t))

dt[, trait := factor(trait, levels = traits)]

means <- dt[, .(mu = mean(log10(ratio))), trait]



ggplot(dt, aes(x = log10(ratio), fill = trait)) +
    geom_density() +
    geom_vline(data = means,
               aes(xintercept = mu),
               linetype = "dashed") +
    labs(x = "Ratio",
         y = "Density") +
    facet_wrap(~ trait,
               nrow = 1,
               labeller = labeller(
                   trait = c(sus = "Susceptibility",
                             inf = "Infectivity",
                             tol = "Tolerance"))) +
    theme_bw() +
    theme(legend.position = "off")
