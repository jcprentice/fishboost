{
    library(data.table)
    library(stringr)
    library(purrr)
}

# Make a table with covariances with mean, ESS, and GR

# dataset <- "fb-final"; scens <- 1:8
dataset <- "fb-test"; scens <- 1:8


x <- map(scens, \(i) {
    # i <- 1
    res_file <- str_glue("datasets/{dataset}/results/scen-{i}-1.rds")
    if (!file.exists(res_file)) return(x)
    
    res <- readRDS(res_file)
    PEs <- res$parameter_estimates
    params <- res$params
    
    scenario <- params$scenario
    trial <- str_split_i(params$setup, "_", 2)
    
    PEs[!str_starts(parameter, "Group effect") & parameter != "latent_period",
              .(parameter, mean, ESS = as.numeric(ESS), GR, scenario, trial,
                trial_fe = params$trial_fe,
                donor_fe = params$donor_fe,
                txd_fe = params$txd_fe,
                weight_fe = params$weight_fe)]
}) |>
    rbindlist() |>
    setorder(parameter, scenario)

# x[parameter == "cov_G_ss", parameter := "cov_G (sus)"]
# x[parameter == "cov_G_ii", parameter := "cov_G (inf)"]
# x[parameter == "cov_G_tt", parameter := "cov_G (mort)"]

# What values do we get with the best convergence?
x[, mean := round(mean, 3)]
x[str_starts(parameter, "cov_G")]
x[str_starts(parameter, "cov_G")][order(GR), .SD, parameter]
x[str_starts(parameter, "cov_G")][order(ESS), .SD, parameter]

# Min and max covariances
y <- x[str_starts(parameter, "cov_G"),
  .(min = round(min(mean), 2), max = round(max(mean), 2)),
  .(parameter, trial)] 
y[, range := str_c("(", min, ", ", max, ")")]
y[, c("min", "max") := NULL]
y
