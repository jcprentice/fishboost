{
    library(data.table)
    library(stringr)
    library(purrr)
}

calc_fixed_effects <- function(dataset = "fb-final",
                               scen = 1,
                               tidy = TRUE) {
    
    # dataset <- "fb-final"; scen <- 1
    
    f <- str_glue("datasets/{dataset}/results/scen-{scen}-1.rds")
    PEs <- readRDS(f)$parameter_estimates
    
    x <- PEs[str_detect(parameter, "beta|period|trial|donor|txd|weight"), .(parameter, mean)]
    FEs <- setNames(x$mean, x$parameter)
    rm(x)
    
    # We need this to ensure donors have been properly allocated
    f <- str_glue("datasets/{dataset}/data/scen-{scen}-1-data.tsv")
    data <- fread(f)
    
    cols <- c("group", "trial_fe", "donor_fe", "txd_fe", "weight_fe", "weight1_fe", "weight2_fe")
    cols <- intersect(cols, names(data))
    
    data <- data[sdp == "progeny", ..cols]
    
    ptrial <- data[, mean(trial_fe)]
    pdonor <- data[, mean(donor_fe)]
    ptxd   <- data[, mean(txd_fe)]
    
    # weights should have mean 0, but we need the cutoffs
    wt_cuts <- if ("weight1_fe" %in% cols) {
        data[, .(W1_1 = quantile(weight1_fe, 1/3),
                 W1_2 = quantile(weight1_fe, 2/3),
                 W2_1 = quantile(weight2_fe, 1/3),
                 W2_2 = quantile(weight2_fe, 2/3)), trial_fe]
    } else if ("weight_fe" %in% cols) {
        data[, .(W_1 = quantile(weight_fe, 1/3),
                 W_2 = quantile(weight_fe, 2/3)), trial_fe]
    } else {
        0
    }
    
    rm(data)
    
    vals <- data.table(Category = "",
                       status = rep(c("Donor", "Recipient"), 2),
                       trial = rep(1:2, each = 2),
                       inf_FE = 0, beta = 0, lat_FE = 0, LP = 0, det_FE = 0, DP = 0, tol_FE = 0, RP = 0)
    
    walk(1:4, \(i) {
        walk(c("inf_FE", "lat_FE", "det_FE", "tol_FE"), \(trait) {
            t1   <- str_sub(trait, 1, 1)
            donor_x <- str_c("donor_", t1)
            trial_x <- str_c("trial_", t1)
            txd_x   <- str_c("txd_", t1)
            vals[i, (trait) := ((status == "Donor") - pdonor) * FEs[[donor_x]] + 
                     ((trial == 2) - ptrial) * FEs[[trial_x]] +
                     ((trial == 2 & status == "Donor") - ptxd) * FEs[[txd_x]]]
            
            if ("weight" %in% cols) {
                weight_x <- str_c("weight_", t1)
                vals[i, (trait) := get(trait) + FEs[[weight_x]]]
            }
        })
    })
    
    vals[, `:=`(Category = str_c(status, " ", trial),
                beta = exp(inf_FE) * FEs[["beta"]],
                LP   = exp(lat_FE) * FEs[["latent_period"]],
                DP   = exp(det_FE) * FEs[["detection_period"]],
                RP   = exp(rec_FE) * FEs[["removal_period"]])]
    vals[, c("status", "trial") := NULL]
 
    if (tidy) {
        vals[, `:=`(inf_FE = round(inf_FE, 1), beta = round(beta, 3),
                    lat_FE = round(lat_FE, 1), LP = round(LP, 1),
                    det_FE = round(det_FE, 1), DP = round(DP, 1),
                    tol_FE = round(tol_FE, 1), RP = round(RP, 1))]
        setnames(vals, c("LP", "DP", "RP"),
                 c("LP (days)", "DP (days)", "RP (days)"))
    }
    
    vals
}

# dataset <- "fb-parasites4"; scen <- 3
# dataset <- "fb-fes2"; scen <- 3
dataset <- "fb-fes3"; scen <- 3

vals <- calc_fixed_effects(dataset, scen)
message(str_glue("dataset '{dataset} / s{scen}'"))
print(vals)
