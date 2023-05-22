library(data.table)
library(glue)

data_set <- "fb-parasites2"

parasites <- c("S", "SE", "SEI", "unused")
trial <- c(1, 2, 12)
num_scens <- length(parasites) * length(trial)

# Expected periods for each category
X <- data.table(scenario = 1:num_scens,
                parasites = rep(parasites,
                                each = length(trial)),
                trial = rep(trial, length(parasites)))

pars <- c("Latent_Period", "Detection_Period",
          "Donor_L", "Trial_L", "TxD_L", "Donor_D", "Trial_D", "TxD_D")

for (i in 1:num_scens) {
    load(glue("results/{data_set}/scen-{i}-1.RData"))
    for (Par in pars) {
        par <- tolower(Par);
        if (par %in% parameter_estimates$parameter) {
            pe <- parameter_estimates[parameter == par, round(mean, 1)]
            X[scenario == i, (Par) := pe]
        }
    }
}

setnames(X, pars[1:2], c("LP", "DP"))

for (scen in 1:num_scens) {
    # Load FEs
    data <- fread(glue("data/{data_set}/scen-{scen}-1-data.tsv"))[sdp == "progeny"]
    
    has_trial <- "trial_fe" %in% names(data)
    
    if (has_trial) {
        # Turn into matrix
        ml <- md <- as.matrix(data[, .(trial_fe, donor_fe, txd_fe)])
        
        # Get idx of each unique type of individual
        idxs <- which(!duplicated(ml))
        ntypes <- length(idxs)
        
        # Convert into posterior values
        ml <- ml %*% diag(as.numeric(X[scenario == scen, .(Trial_L, Donor_L, TxD_L)]), nrow = 3, ncol = 3)
        md <- md %*% diag(as.numeric(X[scenario == scen, .(Trial_D, Donor_D, TxD_D)]), nrow = 3, ncol = 3)
        
        # Subtract means
        ml <- ml - matrix(1, nrow(ml), ncol(ml)) %*% diag(colMeans(ml), nrow = 3, ncol = 3)
        md <- md - matrix(1, nrow(md), ncol(md)) %*% diag(colMeans(md), nrow = 3, ncol = 3)
        
        # Extract typical unique individual and add expected values to X
        rls <- rowSums(ml)[idxs]
        rds <- rowSums(md)[idxs]
        LPs <- c(exp(rls) * X[scenario == scen, LP],
                 exp(rds) * X[scenario == scen, DP])
        
        names(LPs) <- glue("{DR}{trial}_{period}",
                           DR = rep(c("Donor", "Recip"), 4),
                           trial = rep(1:2, each = 2, 2),
                           period = rep(c("LP", "DR"), each = 4))
    } else {
        # Turn into matrix
        ml <- md <- as.matrix(data[, .(donor_fe)])
        
        # Get idx of each unique type of individual
        idxs <- which(!duplicated(ml))
        ntypes <- length(idxs)
        
        # Convert into posterior values
        ml <- ml * X[scenario == scen, Donor_L]
        md <- md * X[scenario == scen, Donor_D]
        
        # Subtract means
        ml <- ml - colMeans(ml)
        md <- md - colMeans(md)
        
        # Extract typical unique individual and add expected values to X
        rls <- rowSums(ml)[idxs]
        rds <- rowSums(md)[idxs]
        LPs <- c(exp(rls) * X[scenario == scen, LP],
                 exp(rds) * X[scenario == scen, DP])
        
        names(LPs) <- glue("{DR}{trial}_{period}",
                           DR = rep(c("Donor", "Recip"), 2),
                           trial = data$trial[[1]],
                           period = rep(c("LP", "DP"), each = 2))
    }
    
    for (lp in names(LPs)) {
        X[scenario == scen, (lp) := round(LPs[[lp]], 3)]
    }
}

message("Data set: ", data_set)
message("Latent Periods")
print(X[, .(scenario, parasites, trial, LP,
            Donor_L, Trial_L, TxD_L,
            Donor1_LP, Recip1_LP, Donor2_LP, Recip2_LP)])
message("Detection Periods")
print(X[, .(scenario, parasites, trial, DP,
            Donor_D, Trial_D, TxD_D,
            Donor1_DP, Recip1_DP, Donor2_DP, Recip2_DP)])
