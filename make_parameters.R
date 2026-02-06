{
    library(data.table)
    library(stringr)
    
    source("widen_priors.R")
    source("utils.R")
}

#' Generate a list of parameters
#'
#' @param model_type Epidemic model: "SIR", "SIDR", "SEIR", "SEIDR"
#' @param name Filename for data
#' @param dataset Folder suffix for data (default "")
#' @param scenario Integer for scenario
#' @param replicate Integer for replicate
#' @param setup Population layout
#' @param use_traits Which traits to use (this is going to be "clever") use
#'   "all" or "none", or a subset using the first letter of each trait
#'   "sildt" = "sus", "inf", "lat", "det", "tol"
#' @param vars Variances: diagonal of cov matrix, either
#'   - a single value, e.g. 1.0
#'   - a list wrapped vector of compatible length, e.g. list(1, 1.5, 0, 0, 0.5)
#'   - a list wrapped named vector e.g. list(sus = 1, inf = 1.5, tol = 0.5, default = 0)
#' @param cors Correlations for Sigma matrix, either
#'   - a single numerical value, e.g. 0.2
#'   - a named list-wrapped vector, e.g. list(si = 0.3, st = -0.1, it = 0.2, default = 0)
#'   - a named list of compatible with the lower triangle order, e.g.
#'     list(si, sl, sd, st, il, id, it, ld, lt, dt)
#' @param group_layout How to arrange individuals into groups: "random",
#'   "striped", "family", "fishboost"
#' @param group_effect Include group effect if >= 0, simulate if appropriate and
#'   check, ignore if < 0
#' @param trial_fe List of fixed effects to check for donors, trial, and weight.
#'   Should also simulate these if not using FB dataset. A string with the first
#'   letters of traits, e.g. "lidt"
#' @param donor_fe
#' @param txd_fe
#' @param weight_fe
#' @param weight_is_nested should weight be nested across trials? (default TRUE),
#' @param sim_new_data Simulate new data in "r" or "bici". Alternatively, "no"
#'   for using FB data, or "etc_sim" or "etc_inf" for getting data from an
#'   existing dataset (specified in patches)
#' @returns a list containing all the parameters

make_parameters <- function(
        model_type = "SEIDR", name = "scen-1-1", dataset = "testing",
        scenario = 1L, replicate = 1L, setup = "fb_12_rpw", use_traits = "sit",
        vars = 1.0, cors = 0.2, group_layout = "fishboost", group_effect = -1,
        trial_fe = "none", donor_fe = "none", txd_fe = "none",
        weight_fe = "none", weight_is_nested = TRUE, sim_new_data = "bici"
) {
    message("Setting parameters ...")
    
    # Fix input NAs ----
    
    # Handy in case of debugging
    if (FALSE) {
        model_type <- "SEIDR"; name <- "scen-1-1"; dataset <- "testing";
        scenario <- 1L; replicate = 1L; setup <- "fb_12_rpw"; use_traits <- "sit";
        vars <- 1; cors <- 0.2; group_layout <- "fishboost"; group_effect <- -1;
        trial_fe <- "ildt"; donor_fe <- "ildt"; txd_fe <- "ildt";
        weight_fe <- "sildt"; weight_is_nested <- TRUE; sim_new_data <- "bici";
    }
    
    if (FALSE) {
        model_type <- protocol$model_type; dataset <- protocol$dataset; name <- protocol$name;
        scenario <- protocol$scenario; replicate = protocol$scenario; use_traits <- protocol$use_traits;
        vars <- protocol$vars; cors <- protocol$cors; setup <- protocol$setup;
        group_layout <- protocol$group_layout; group_effect <- protocol$group_effect;
        trial_fe <- protocol$trial_fe; donor_fe <- protocol$donor_fe; txd_fe <- protocol$txd_fe;
        weight_fe <- protocol$weight_fe; weight_is_nested <- protocol$weight_is_nested;
        sim_new_data <- protocol$sim_new_data;
    }
    
    # Fix any parameters incorrectly assigned as NULL
    {
        model_type    <- model_type %||% "SEIDR"
        name          <- name %||% "scen-1-1"
        dataset       <- dataset %||% "testing"
        scenario      <- scenario %||% 1L
        replicate     <- replicate %||% 1L
        setup         <- setup %||% "fb_12_rpw"
        use_traits    <- use_traits %||% "sit"
        vars          <- vars %||% 1
        cors          <- cors %||% 0.2
        group_layout  <- group_layout %||% "fishboost"
        trial_fe      <- trial_fe %||% "none"
        donor_fe      <- donor_fe %||% "none"
        txd_fe        <- txd_fe %||% "none"
        weight_fe     <- weight_fe %||% "none"
        sim_new_data  <- sim_new_data %||% "bici"
    }
    
    # Description
    description <- "Basic test model"

    label <- str_c("s", scenario)
    
    # Set output directories ----
    
    # Directories should be something like "dataset/data", "dataset/results", or
    # "testing" if just testing (and convert to basic string)
    if (dataset == "") {
        dataset <- "testing"
    }
    
    {
        base_dir    <- c(str_glue("datasets/{dataset}"))
        gfx_dir     <- c(str_glue("{base_dir}/gfx"))
        data_dir    <- c(str_glue("{base_dir}/data"))
        results_dir <- c(str_glue("{base_dir}/results"))
        meta_dir    <- c(str_glue("{base_dir}/meta"))
        output_dir  <- c(str_glue("{data_dir}/{name}-out"))
        states_dir  <- c(str_glue("{output_dir}/states"))
        config      <- c(str_glue("{data_dir}/{name}"))
    }
    
    # Population setup ----
    message(str_glue("- Setup is: '{setup}'\n",
                     "- Group layout is: '{group_layout}'"))
    
    # Note: the FB dataset isn't balanced
    # otherwise dpsire = dams per sire, ppdam = progeny per dam
    tmp <- switch(
        setup,           # sires dams  progeny groups trials I0
        "fb_12"        = c(29,   25,   1775,   71,    2,     5),
        "fb_1"         = c(14,   14,   875,    35,    1,     5),
        "fb_2"         = c(18,   14,   900,    36,    1,     5),
        "fb_12_drop71" = c(29,   25,   1750,   70,    2,     5),
        "fb_1_drop71"  = c(14,   14,   875,    35,    1,     5),
        "fb_2_drop71"  = c(18,   14,   875,    35,    1,     5),
        "fb_12_rpw"    = c(28,   25,   1750,   70,    2,     5),
        "fb_1_rpw"     = c(14,   14,   875,    35,    1,     5),
        "fb_2_rpw"     = c(17,   14,   875,    35,    1,     5),
        "chris"        = c(100,  2000, 2000,   200,   1,     5),
        "small"        = c(3,    6,    12,     4,     1,     1),
        "single"       = c(10,   20,   500,    1,     1,     5),
        "multiple"     = c(10,   20,   500,    20,    1,     5),
                         c(28,   25,   1750,   70,    2,     5)) |>
        as.integer() |> 
        setNames(c("nsires", "ndams", "nprogeny", "ngroups", "ntrials", "I0")) |>
        as.list()
    
    nsires <- tmp$nsires; ndams <- tmp$ndams; nprogeny <- tmp$nprogeny
    ngroups <- tmp$ngroups; ntrials <- tmp$ntrials; I0 <- tmp$I0
    
    
    iceil <- function(x) as.integer(ceiling(x))
    
    # Derived numbers
    {
        nparents <- nsires + ndams
        ntotal <- nprogeny + nparents
        dpsire <- iceil(ndams / nsires)
        ppdam <- iceil(nprogeny / ndams)
        group_size <- iceil(nprogeny / ngroups)
    }
    
    # No trial FEs if only one trial
    if (ntrials == 1) {
        trial_fe <- txd_fe <- sim_trial_fe <- sim_txd_fe <- ""
        weight_is_nested <- FALSE
    }
    
    # Traits ----
    message(str_glue("- Model type is: '{model_type}'"))
    
    # Compartments are just the letters in model_type (ignore repeated S)
    compartments <- uniq_chars(model_type)
    
    all_traits <- c(s = "sus", i = "inf", l = "lat", d = "det", t = "tol")
    
    # Set up which traits a model will use, if traits are +vely or -vely
    # correlated, what the event times are called, and which parameters we need
    # to make SIRE aware of.
    switch(model_type,
           "SEIDR" = {
               model_traits <- all_traits
               timings <- c("Tinf", "Tinc", "Tsym", "Tdeath")
           }, "SIDR" = {
               model_traits <- all_traits[c("s", "i", "d", "t")]
               timings <- c("Tinf", "Tsym", "Tdeath")
           }, "SEIR" = {
               model_traits <- all_traits[c("s", "i", "l", "t")]
               timings <- c("Tinf", "Tsym", "Tdeath")
           }, "SIR" = {
               model_traits <- all_traits[c("s", "i", "t")]
               timings <- c("Tinf", "Tdeath")
           }, "SIS" = {
               model_traits <- all_traits[c("s", "i", "t")]
               timings <- c("Tinf", "Tdeath")
           }, "SI" = {
               model_traits <- all_traits[c("s", "i")]
               timings <- c("Tinf")
           }, {
               stop("- Unknown model!")
        })
    
    ## Genetic and Environmental covariances ----
    
    # Likely using only a subset of traits
    use_traits <- if (use_traits %in% c("", "none", NA)) {
        ""
    } else {
        str_to_lower(use_traits)
    }
    
    n_traits <- length(model_traits)
    n_traits_used <- str_length(use_traits)
    traits_used <- model_traits[str_chars(use_traits)]
    
    if (n_traits_used > 0) {
        message("- ", n_traits_used, " GE traits: ", str_flatten_comma(traits_used))
    } else {
        message("- 0 GE traits used")
    }
    
    ## Covariance matrices ----
    
    out <- make_matrices(all_traits, use_traits, vars, cors)
    Sigma_G <- Sigma_E <- out$Sigma
    cov_G <- cov_E <- out$cov
    
    # This should expand vars and cors to be the full matrix values, possibly
    # with replicates
    vars <- Sigma_G |> diag()
    cors <- Sigma_G[lower.tri(Sigma_G)] |> setNames(out$cor_names)
    
    ge_opts <- "" # c("gt_only", "pt_only", "e1", "no_ev_xyz")
    
    
    # Epidemic parameters ----
    
    ## Infection coefficient ----
    r_beta <- 0.5
    
    # Infectivity model
    # 1: standard
    # 2: I is 10% as infectious as D
    # 3: donors are 10% as infectious as recipients
    # 4: donors are `id_ratio` as infectious as recipients
    inf_model <- 1L
    
    # Infectivity ratio
    inf_ratio <- if (inf_model == 1) 1 else 0.1
    
    # Note for Gamma dist:
    # shape = mean^2 / var
    # rate  = mean / var
    
    ## Compartment Periods ----
    
    ### Latency period E->I ----
    
    # Calculated for an SEIDR model based on donors in trial 1 (+ 0.5)
    # fitdist(Tsym + 0.5, "gamma")
    latent_period <- 10
    LP_shape <- 1 # 1.46 # 2.0 # 10.0
    LP_scale <- latent_period / LP_shape # 0.232 # LP_shape / latent_period
    LP_dist <- "exp"
    
    
    ### Detection period I->D ----
    
    # This is currently identical to the latency
    detection_period <- 10
    DP_shape <- 1 # 1.46 # 2.0 # 10.0
    DP_scale <- detection_period / DP_shape # 0.232 # LP_shape / latent_period
    DP_dist <- "exp"
    
    
    ### Removal ----
    
    # fitdist(RP + 0.5, "gamma")
    removal_period <- 10
    RP_shape <- 1
    RP_scale <- removal_period / RP_shape # 0.149
    RP_dist <- "exp"
    
    
    ## Fixed effects ----
    
    # prevent fread interpreting these as NA
    if (trial_fe  %in% c("none", "", NA)) trial_fe  <- ""
    if (donor_fe  %in% c("none", "", NA)) donor_fe  <- ""
    if (txd_fe    %in% c("none", "", NA)) txd_fe    <- ""
    if (weight_fe %in% c("none", "", NA)) weight_fe <- ""
    
    # when simulating, copy these values, or override them after params is
    # created if they should be different
    sim_trial_fe  <- trial_fe
    sim_donor_fe  <- donor_fe
    sim_txd_fe    <- txd_fe
    sim_weight_fe <- weight_fe
    
    # A default set of FEs
    
    #                    sus   inf   lat   det   tol
    fe_vals <- matrix(c(+0.0, -0.8, -2.9, +2.9, -0.2,  # trial
                        +0.0, -2.1, -1.2, +1.2, -1.4,  # donor
                        +0.0, -2.4, +3.3, +1.4, +1.7,  # txd
                        -0.1, +0.2, -0.8, +0.7, +0.2,  # weight
                        +0.1, -0.6, -0.2, +1.6, +0.3,  # weight1
                        -0.2, +1.6, +4.0, +0.1, +0.1), # weight2
                      ncol = 5, nrow = 6, byrow = TRUE,
                      dimnames = list(fe = c("trial", "donor", "txd", "weight", "weight1", "weight2"),
                                      traits = unname(all_traits)))
    
    fe_vals[] <- 0
    
    ## Map traits and FEs ----
    
    # "xxxxx" will be copied onto "sildt", so "sittt" means that lat, det, tol
    # will all use the same tol value. sim_ is for simulating new data in R,
    # otherwise it's what SIRE or BICI is passed
    {
        sim_link_traits <- "sildt"
        sim_link_trial  <- "sildt"
        sim_link_donor  <- "sildt"
        sim_link_txd    <- "sildt"
        sim_link_weight <- "sildt"
        
        link_traits <- "sildt"
        link_trial  <- "sildt"
        link_donor  <- "sildt"
        link_txd    <- "sildt"
        link_weight <- "sildt"
        
        # Shape parameters for lat, det, tol
        sim_link_shapes <- "ldt"
        link_shapes <- "ldt"
    }    
    
    ## R0 ----
    
    # This should target R0 around 5?
    R0 <- r_beta * (removal_period + if ("D" %in% compartments) detection_period else 0)
    
    
    # MCMC settings ----
    
    bici_cmd <- "inf"
    
    nchains <- 16L
    
    ## Samples ----
    nsample  <- 1e4 # total no. of samples SIRE should take
    burnprop <- 0.2 # burn-in proportion of chain 
    thinto   <- 1e4 # no. of samples to output
    burnin   <- nsample * burnprop
    thin     <- max(nsample / thinto, 1)
    
    # This should be tweaked to give ~ 25% anneal time
    nsample_per_gen <- max(3e-3 * nsample, 1)
    
    sample_states <- 0L # how many states (per chain) to sample
    ie_output <- "true" # if ss > 0, include IE samples in extended_trace_combine.tsv
    
    # Should SIRE use the PAS method (N chains for a single call) or regular MCMC?
    algorithm    <- "pas" # "mcmc"
    anneal       <- "on"  # not used in "pas"
    anneal_power <- 4     # not used in "pas"
    
    phi <- 1.0
    
    # Priors ----
    
    ge <- max(group_effect, 0)

    # Do we use a Jeffreys or a uniform prior for covariance
    cov_prior <- c("jeffreys", "uniform")[[1]]
    single_prior <- c("inverse", "uniform")[[1]]
    
    
    priors <- rowwiseDT(
        parameter=,         type=,        val1=, val2=, true_val=,
        
        "beta",             single_prior, 0,     2.0,   r_beta,
        "latent_period",    single_prior, 0,     40,    latent_period,
        "detection_period", single_prior, 0,     160,   detection_period,
        "removal_period",   single_prior, 0,     20,    removal_period,
        "LP_shape",         "uniform",    0.5,   5,     LP_shape,
        "DP_shape",         "uniform",    0.5,   5,     DP_shape,
        "RP_shape",         "uniform",    0.5,   5,     RP_shape,
        
        # parameter                 type          val1 val2 true_val
        "beta_Tr1",                 single_prior, 0,   2,   r_beta,
        "beta_Tr2",                 single_prior, 0,   2,   r_beta,
        "latent_period_Tr1,Don",    single_prior, 0,   50,  latent_period,
        "latent_period_Tr1,Rec",    single_prior, 0,   50,  latent_period,
        "latent_period_Tr2,Don",    single_prior, 0,   120, latent_period,
        "latent_period_Tr2,Rec",    single_prior, 0,   50,  latent_period,
        "detection_period_Tr1,Don", single_prior, 0,   50,  detection_period,
        "detection_period_Tr1,Rec", single_prior, 0,   50,  detection_period,
        "detection_period_Tr2,Don", single_prior, 0,   160, detection_period,
        "detection_period_Tr2,Rec", single_prior, 0,   50,  detection_period,
        "removal_period_Tr1,Don",   single_prior, 0,   30,  removal_period,
        "removal_period_Tr1,Rec",   single_prior, 0,   30,  removal_period,
        "removal_period_Tr2,Don",   single_prior, 0,   30,  removal_period,
        "removal_period_Tr2,Rec",   single_prior, 0,   30,  removal_period,
        "infrat",                   "uniform",    0,    1,  inf_ratio,

        # parameter type       val1 val2    true_val
        "sigma",    "inverse", 0,   2 * ge, ge,
        
        # parameter type       val1  val2 true_val
        "cov_G_ss", "uniform",    0,   4, Sigma_G["sus", "sus"],
        "cov_G_ii", "uniform",    0,   4, Sigma_G["inf", "inf"],
        "cov_G_ll", "uniform",    0,   2, Sigma_G["lat", "lat"],
        "cov_G_dd", "uniform",    0,   2, Sigma_G["det", "det"],
        "cov_G_tt", "uniform",    0, 1.5, Sigma_G["tol", "tol"],
        "r_G_si",   "uniform", -0.9, 0.9, Sigma_G["sus", "inf"],
        "r_G_sl",   "uniform", -0.9, 0.9, Sigma_G["sus", "lat"],
        "r_G_sd",   "uniform", -0.9, 0.9, Sigma_G["sus", "det"],
        "r_G_st",   "uniform", -0.9, 0.9, Sigma_G["sus", "tol"],
        "r_G_il",   "uniform", -0.9, 0.9, Sigma_G["inf", "lat"],
        "r_G_id",   "uniform", -0.9, 0.9, Sigma_G["inf", "det"],
        "r_G_it",   "uniform", -0.9, 0.9, Sigma_G["inf", "tol"],
        "r_G_ld",   "uniform", -0.9, 0.9, Sigma_G["lat", "tol"],
        "r_G_lt",   "uniform", -0.9, 0.9, Sigma_G["lat", "tol"],
        "r_G_dt",   "uniform", -0.9, 0.9, Sigma_G["det", "tol"],
        "cov_E_ss", "uniform",    0,   4, Sigma_E["sus", "sus"],
        "cov_E_ii", "uniform",    0,   4, Sigma_E["inf", "inf"],
        "cov_E_ll", "uniform",    0,   4, Sigma_E["lat", "lat"],
        "cov_E_dd", "uniform",    0,   4, Sigma_E["det", "det"],
        "cov_E_tt", "uniform",    0,   1, Sigma_E["tol", "tol"],
        "r_E_si",   "uniform", -0.9, 0.9, Sigma_E["sus", "inf"],
        "r_E_sl",   "uniform", -0.9, 0.9, Sigma_E["sus", "lat"],
        "r_E_sd",   "uniform", -0.9, 0.9, Sigma_E["sus", "det"],
        "r_E_st",   "uniform", -0.9, 0.9, Sigma_E["sus", "tol"],
        "r_E_il",   "uniform", -0.9, 0.9, Sigma_E["inf", "lat"],
        "r_E_id",   "uniform", -0.9, 0.9, Sigma_E["inf", "det"],
        "r_E_it",   "uniform", -0.9, 0.9, Sigma_E["inf", "tol"],
        "r_E_ld",   "uniform", -0.9, 0.9, Sigma_E["lat", "tol"],
        "r_E_lt",   "uniform", -0.9, 0.9, Sigma_E["lat", "tol"],
        "r_E_dt",   "uniform", -0.9, 0.9, Sigma_E["det", "tol"],
        
        # parameter  type       val1 val2 true_val
        "trial_i",   "uniform", -5, 3, fe_vals["trial", "inf"],
        "trial_l",   "uniform", -5, 3, fe_vals["trial", "lat"],
        "trial_d",   "uniform", -4, 4, fe_vals["trial", "det"],
        "trial_t",   "uniform", -4, 4, fe_vals["trial", "tol"],
        
        "donor_i",   "uniform", -7, 1, fe_vals["donor", "inf"],
        "donor_l",   "uniform", -5, 3, fe_vals["donor", "lat"],
        "donor_d",   "uniform", -4, 4, fe_vals["donor", "det"],
        "donor_t",   "uniform", -5, 3, fe_vals["donor", "tol"],
        
        "txd_i",     "uniform", -5, 3, fe_vals["txd", "inf"],
        "txd_l",     "uniform", -1, 7, fe_vals["txd", "lat"],
        "txd_d",     "uniform", -2, 6, fe_vals["txd", "det"],
        "txd_t",     "uniform", -2, 6, fe_vals["txd", "tol"],
        
        "weight_s",  "uniform", -4, 4, fe_vals["weight", "sus"],
        "weight_i",  "uniform", -4, 4, fe_vals["weight", "inf"],
        "weight_l",  "uniform", -4, 4, fe_vals["weight", "lat"],
        "weight_d",  "uniform", -4, 4, fe_vals["weight", "det"],
        "weight_t",  "uniform", -4, 4, fe_vals["weight", "tol"],
        
        "weight1_s", "uniform", -3, 3, fe_vals["weight1", "sus"],
        "weight1_i", "uniform", -3, 3, fe_vals["weight1", "inf"],
        "weight1_l", "uniform", -3, 3, fe_vals["weight1", "lat"],
        "weight1_d", "uniform", -3, 3, fe_vals["weight1", "det"],
        "weight1_t", "uniform", -3, 3, fe_vals["weight1", "tol"],
        
        "weight2_s", "uniform", -3, 3, fe_vals["weight2", "sus"],
        "weight2_i", "uniform", -1, 7, fe_vals["weight2", "inf"],
        "weight2_l", "uniform",  0, 8, fe_vals["weight2", "lat"],
        "weight2_d", "uniform", -3, 3, fe_vals["weight2", "det"],
        "weight2_t", "uniform", -3, 3, fe_vals["weight2", "tol"]
    )
    
    # Roughly fix some priors
    widen_priors(priors)
    
    priors[, use := FALSE]
    
    # This will be all be redone with set_use_flags(), but make_parameters()
    # should give as correct a result as possible.
    
    
    used <- str_1st(model_traits)
    
    use_parameters <- c(
        "beta",
        if ("l" %in% used) "latent_period",
        if ("d" %in% used) "detection_period",
        if ("t" %in% used) "removal_period",
        if ("l" %in% used && LP_dist == "gamma") "LP_shape",
        if ("d" %in% used && DP_dist == "gamma") "DP_shape",
        if ("p" %in% used && RP_dist == "gamma") "RP_shape",
        if (group_effect > 0) "sigma"
    )
    priors[parameter %in% use_parameters, use := TRUE]
    
    GE_traits <- used |>
        intersect(str_chars(use_traits)) |>
        intersect(str_chars(link_traits))
    
    pwalk(expand.grid(x = used, y = used), \(x, y) {
        priors[str_ends(parameter, str_glue("_[GE]_{x}{y}")),
               use := all(c(x, y) %in% GE_traits)]
    })
    
    fe_str <- function(base) {
        # fe_str("trial") -> "trial_[ildt]"
        # fe_str("weight") -> "weight[12]_[ildt]"
        
        base_str <- get(str_glue("{base}_fe"))
        link_str <- get(str_glue("link_{base}"))
        
        if (base_str %in% c("", "none")) return("XYZ")
        
        used |>
            intersect(str_chars(base_str)) |>
            intersect(str_chars(link_str)) |>
            str_flatten() |>
            str_glue("{base}{wt}_[{fe}]",
                     fe = _,
                     wt = if (base == "weight" && weight_is_nested)
                         "[12]" else "")
    }
    
    priors[str_detect(parameter, fe_str("trial")),  use := TRUE]
    priors[str_detect(parameter, fe_str("donor")),  use := TRUE]
    priors[str_detect(parameter, fe_str("txd")),    use := TRUE]
    priors[str_detect(parameter, fe_str("weight")), use := TRUE]
    
    
    # Additional settings ----
    
    popn_format <- c("intervals", "times")[[2]]

    seed <- 0

    # Number of times to simulate in BICI
    nreps <- 50
    
    # What we pass to SIRE 2.0, choice of
    # "Tinf": actual infection time
    # "Tsym": time of symptoms
    # "estimated_Tinf_from_donors": Tsym - mean incubation period of donors
    # "estimated_Tinf_per_individual": Tsym with all gaps covered
    pass_Tsym <- "Tsym"
    
    # pass event times, pass the last N of timings (Tinf, Tinc, Tsym, Tdeath)
    pass_events <- c("Tsym", "Tdeath") # pass time Tsym and Tdeath

    # time step between data
    time_step <- 0
    # time step for inference
    time_step_bici <- 1
    
    censor <- 1.0
    
    ## Fixes for FB data set
    fix_donors <- c() # c("time", "no_Tsym_survivors", "set_to_R")
    # All donors who fail to show symptoms by this point are demoted
    t_demote <- c(20, 80)
    # if Tsym == Tdeath, then bump Tsym back by 1/2
    fix_eq_time <- TRUE

    
    # Patch data with values from dataset / scenario posterior
    patch_dataset <- ""
    patch_name <- "scen-1-1"
    # Use mean / median / randomly sample from posterior
    patch_type <- c("mean", "median", "sampled")[[1]]
    # Use states to patch EBVs and parameters (otherwise use trace files)
    patch_state <- FALSE
    
    # Use weight if available
    use_weight <- "log" # "log" or "linear"
    
    # Use inverse H matrix, GRM, or none
    use_grm <- "" # {A,H} + {,_inv} + {,_nz}
    
    # How to build popn?
    # "pedigree" / "grm" = generate new values based on relationship
    # "posterior" = patch from posteriors given in "patch_xyz"
    traits_source <- c("pedigree", "grm", "posterior")[[1]]
    
    # Show plots during run
    show_plots <- TRUE
    
    # show messages
    msgs <- TRUE
    
    # Show details
    DEBUG <- FALSE
    
    
    # Create params list ----
    #
    # This will be a list for passing to functions
    params <- mget(
        # Inputs
        c("description", "dataset", "scenario", "name", "label", "replicate", "seed",
          "data_dir", "results_dir", "gfx_dir", "meta_dir", "output_dir", "states_dir", "config",
          "model_type", "sim_new_data", "setup", "vars", "cors", "group_layout", "use_traits",
          # Population parameters
          "ntrials", "nsires", "ndams", "nparents", "nprogeny", "ngroups", "ntotal",
          "dpsire", "ppdam", "group_size", "I0",
          # Model traits
          "compartments", "timings", "all_traits", "model_traits",
          "n_traits", "Sigma_E", "Sigma_G", "cov_G", "cov_E", "ge_opts",
          "sim_link_traits", "sim_link_trial", "sim_link_donor", "sim_link_txd",
          "sim_link_weight", "sim_link_shapes",
          "link_traits", "link_trial", "link_donor", "link_txd", "link_weight", "link_shapes",
          # Main model parameters
          "r_beta", "inf_ratio", "inf_model", "group_effect", "R0", # infection
          "latent_period",    "LP_shape", "LP_scale", "LP_dist", # latency
          "detection_period", "DP_shape", "DP_scale", "DP_dist", # detection
          "removal_period",   "RP_shape", "RP_scale", "RP_dist", # removal
          "fe_vals", "trial_fe", "donor_fe", "txd_fe",
          "sim_trial_fe", "sim_donor_fe", "sim_txd_fe",
          "weight_fe", "sim_weight_fe", "weight_is_nested",
          # MCMC & extra
          "priors", "cov_prior", "single_prior",
          "popn_format", "bici_cmd", "algorithm", "anneal", "anneal_power",
          "nchains", "nsample", "burnprop", "thinto", "nchains",
          "phi", "nsample_per_gen", "sample_states",
          "ie_output", "time_step", "time_step_bici", "censor", "nreps",
          "fix_donors", "t_demote", "fix_eq_time",
          "pass_events", "patch_dataset", "patch_name", "patch_type", "patch_state",
          "use_weight", "use_grm", "traits_source",
          "show_plots", "msgs", "DEBUG"))
    
    params
}


# Handy summary ----
summarise_params <- function(params) {
    with(params, {
        if (sim_new_data == "no") {
            message("Using Fishboost data")
        } else if (sim_new_data == "etc_inf") {
            message("Using data simulated in BICI from prior")
        } else if (sim_new_data == "etc_ps") {
            message("Using data simulated in BICI from posterior")
        } else {
            message(str_glue("Simulating new {model_type} data via {x}",
                             x = str_to_upper(sim_new_data)))
        }

        message(str_glue(
            "- Demography is:",
            "\t{nsires} sires, {ndams} dams, {nprogeny} progeny ({ntotal} total)",
            "\t{ngroups} group{s} (group size {group_size})",
            s = if (ngroups > 1) "s" else "",
            .trim = FALSE, .sep = "\n"
        ))
        
        message(str_glue(
            "- Individual Effects on:",
            str_c("\t", if (use_traits %in% c("", "none")) "none" else
                str_flatten_comma(model_traits[str_chars(use_traits)])),
            .trim = FALSE, .sep = "\n"
        ))
        
        message("- Sigma_G: ")
        capture_message(Sigma_G[model_traits, model_traits])
        
        message(str_glue(
            "- Fixed Effects are:",
            "\tTrial = '{tfe}', Donor = '{dfe}', TxD = '{xfe}', Weight = '{wfe}'",
            tfe = if (trial_fe  %in% c("", "none")) "none" else trial_fe,
            dfe = if (donor_fe  %in% c("", "none")) "none" else donor_fe,
            xfe = if (txd_fe    %in% c("", "none")) "none" else txd_fe,
            wfe = if (weight_fe %in% c("", "none")) "none" else weight_fe,
            .trim = FALSE, .sep = "\n"
        ))
        
        message(str_glue(
            "- Estimated R0:",
            "\t{r0}",
            r0 = signif(R0, 3),
            .trim = FALSE, .sep = "\n"
        ))
        
        message(str_glue(
            "- Running MCMC with:",
            "\t{ns} updates, {th} samples, {burnprop} burnin, {nchains} chains",
            "- BICI script file:",
            "\t'{data_dir}/{name}.bici'",
            "- Results file:",
            "\t'{results_dir}/{name}.rds'",
            ns = format(nsample, scientific = FALSE, big.mark = ","),
            th = format(thinto,  scientific = FALSE, big.mark = ","),
            .trim = FALSE, .sep = "\n"
        ))
    })
}
