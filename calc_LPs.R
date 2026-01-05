library(data.table)
library(purrr)
library(stringr)

dataset <- "fb-final"; scens <- 1:8

pars <- c("latent_period", "detection_period",
          "donor_l", "trial_l", "txd_l", "donor_d", "trial_d", "txd_d")

X <- map_df(scens, \(i) {
    f <- str_glue("datasets/{dataset}/results/scen-{i}-1.rds")
    x <- readRDS(f)$parameter_estimates
    PEs <- x[parameter %in% pars,
             .(parameter, value = round(mean, 1))]
    c(scen = i, setNames(PEs$value, PEs$parameter))
}) |>
    setDT() |>
    setnames(c("latent_period", "detection_period"), c("LP", "DP"))


walk(scens, \(i) {
    # Load FEs
    f <- str_glue("datasets/{dataset}/data/scen-{i}-1-data.tsv")
    data <- fread(f)[sdp == "progeny"]
    
    has_trial <- "trial_fe" %in% names(data)
    
    if (has_trial) {
        # Turn into matrix
        ml <- md <- as.matrix(data[, .(trial_fe, donor_fe, txd_fe)])
        
        # Get idx of each unique type of individual
        idxs <- which(!duplicated(ml))
        ntypes <- length(idxs)
        
        # Convert into posterior values
        ml <- ml %*% diag(as.numeric(X[scen == i, .(trial_l, donor_l, txd_l)]), nrow = 3, ncol = 3)
        md <- md %*% diag(as.numeric(X[scen == i, .(trial_d, donor_d, txd_d)]), nrow = 3, ncol = 3)
        
        # Subtract means
        ml <- ml - matrix(1, nrow(ml), ncol(ml)) %*% diag(colMeans(ml), nrow = 3, ncol = 3)
        md <- md - matrix(1, nrow(md), ncol(md)) %*% diag(colMeans(md), nrow = 3, ncol = 3)
        
        # Extract typical unique individual and add expected values to X
        rls <- rowSums(ml)[idxs]
        rds <- rowSums(md)[idxs]
        LPs <- setNames(c(exp(rls) * X[scen == i, LP],
                          exp(rds) * X[scen == i, DP]) |>
                            round(1),
                        str_glue("{dr}{trial}_{period}",
                                 dr = rep(c("Donor", "Recip"), 4),
                                 trial = rep(1:2, each = 2, 2),
                                 period = rep(c("LP", "DP"), each = 4)))
    } else {
        # Turn into matrix
        ml <- md <- as.matrix(data[, .(donor_fe)])
        
        # Get idx of each unique type of individual
        idxs <- which(!duplicated(ml))
        ntypes <- length(idxs)
        
        # Convert into posterior values
        ml <- ml * X[scen == i, Donor_L]
        md <- md * X[scen == i, Donor_D]
        
        # Subtract means
        ml <- ml - colMeans(ml)
        md <- md - colMeans(md)
        
        # Extract typical unique individual and add expected values to X
        rls <- rowSums(ml)[idxs]
        rds <- rowSums(md)[idxs]
        
        LPs <- setNames(c(exp(rls) * X[scen == i, LP],
                          exp(rds) * X[scen == i, DP]) |>
                            round(1),
                        str_glue("{DR}{trial}_{period}",
                                 DR = rep(c("Donor", "Recip"), 2),
                                 trial = data$trial[[1]],
                                 period = rep(c("LP", "DP"), each = 2)))
    }
    
    walk(names(LPs), ~ X[scenario == i, (.x) := round(LPs[[.x]], 3)])
})

message(str_glue("Data set: '{dataset}'"))
message("Latent Periods")
print(X[, .(scenario, LP,
            donor_l, trial_l, txd_l,
            Donor1_LP, Recip1_LP, Donor2_LP, Recip2_LP)])
message("Detection Periods")
print(X[, .(scenario, DP,
            donor_d, trial_d, txd_d,
            Donor1_DP, Recip1_DP, Donor2_DP, Recip2_DP)])
