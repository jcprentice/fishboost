{
    source("libraries.R")
    source("source_files.R")
    source("crop.R")
    source("set_cors.R")
    library(purrr)
}

## Set parameters ----

cmd_args <- commandArgs(trailingOnly = TRUE)
run_from_script <- length(cmd_args) > 0

if (run_from_script) capture_message(cmd_args)

cv_type <- if (run_from_script) as.numeric(cmd_args[[1L]]) else 2L
cors <- c("none", "none", "sit_pos_cor", "sit_neg_cor", "sit_no_cor")[cv_type]

my_seed <- if (run_from_script) as.numeric(cmd_args[[2L]]) else 1L
set.seed(my_seed)

{
    params <- make_parameters(
        model_type = "SEIDR", # "SIR", "SEIR", "SIDR", or "SEIDR"
        name = "expt1",
        setup = "fb_12", # chris, small, fb_12, fb_1, fb_2, single
        use_traits = "sit", # "all", "none", "sit", "si" etc.
        vars = 1, # c(0.5, 1.2, 0.5)
        cors = cors, # 0.2, "positive", "negative", "mixed", "sit_only" etc.
        group_layout = "fishboost", # "random", "family", "striped", "fishboost"
        trial_fe = "ildt",
        donor_fe = "ildt",
        txd_fe = "ildt",
        weight_fe = "sildt",
        weight_is_nested = TRUE,
        sim_new_data = "r"
    )
    
    params$traits_source <- "pedigree"
    params$patch_dataset <- if (cv_type == 2) "fb-notraits" else "fb-final"
    params$patch_name <- "scen-1-1"
    params$patch_state <- FALSE
    params$patch_type <- "mean"
    params$msgs <- TRUE
    # params$I0 <- 1L
}

# back up parameters
params_base <- copy(params)

params <- params_base |>
    patch_params() |>
    fix_weight_fes() |>
    apply_links()


if (cv_type == 1L) {
    # Do nothing
} else if (cv_type == 2L) {
    # Minimal covariance
    params$Sigma_G[] <- params$Sigma_E[] <- 0
    diag(params$Sigma_G) <- 1e-6
    diag(params$Sigma_E) <- 1e-6
} else {
    # In cases 3-5 keep variances but change covariances
    this_cor <- c(0, 0.2, 0, -0.2, 0)[cv_type]
    params$Sigma_G <- set_cors(params$Sigma_G, this_cor)
    params$Sigma_E <- set_cors(params$Sigma_E, this_cor)
}

popn <- make_pedigree(params) |>
    set_traits(params) |>
    set_groups(params)

# Add in R0
popn[, r0 := sus_g + inf_g + tol_g]

# Set the top percentile to remove
prop <- 0.2

R0s <- map(c(none = "none", sus = "sus_g", r0 = "r0"),
           \(trait) popn |>
               crop(trait, prop) |>
               apply_fixed_effects(params) |>
               simulate_epidemic(params) |>
               get_R0())

# Save R0s
capture_message(unlist(R0s))

if (!dir.exists("r0_selection")) {
        message(" - mkdir r0_selection")
    dir.create("r0_selection")
}

out_file <- str_glue("r0_selection/covar_{cv_type}-{my_seed}.rds")
saveRDS(R0s, file = out_file)

