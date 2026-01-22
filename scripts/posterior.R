{
    library(data.table)
    library(stringr)
    library(purrr)
    library(cowplot)
    source("plotting/plot_posteriors.R")
}

get_posterior <- function(dataset = "fb-final", scen = 1, rep = 1) {
    # dataset <- "fb-test"; scen <- 1; rep <- 1
    # dataset <- "sim-base-inf"; scen <- 1; rep <- 1
    
    {
        base_dir <- str_glue("datasets/{dataset}")
        res_dir  <- str_glue("{base_dir}/results")
        gfx_dir  <- str_glue("{base_dir}/gfx")
        post_dir <- str_glue("{gfx_dir}/posteriors")
        cov_dir  <- str_glue("{gfx_dir}/posteriors_covs")
        
        walk(c(gfx_dir, post_dir, cov_dir), \(d) {
            if (!dir.exists(d)) {
                message(" - mkdir ", d)
                dir.create(d, recursive = TRUE)
            }
        })
    }
    
    plts <- plot_posteriors(dataset = dataset, name = str_glue("{scen}-{rep}"), ci = "hpdi", draw = "density")
    
    
    # Create a plot with subplots aligned by trait
    empty <- ggplot() + theme_classic()
    sildt1 <- c("s", "i", "l", "d", "t")
    sildt2 <- str_c(sildt1, sildt1)
    any_non_empty <- function(x) any(x != "empty")
    
    plts$plts$empty <- empty
    parnames <- names(plts$plts)
    
    cov_pars <- c(str_c("cov_G_", sildt2),
                  "r_G_si", "r_G_st", "empty", "empty", "r_G_it",
                  str_c("cov_E_", sildt2),
                  str_c("cov_P_", sildt2))
    
    # Do we want to include the latent_period here?
    model_pars <- c(
        "sigma",  "beta_Tr1", "LP_Tr1,Don", "DP_Tr1,Don", "RP_Tr1,Don", 
        "infrat", "empty",    "LP_Tr1,Rec", "DP_Tr1,Rec", "RP_Tr1,Rec",
        "sigma",  "beta_Tr2", "LP_Tr2,Don", "DP_Tr2,Don", "RP_Tr2,Don",
        "infrat", "empty",    "LP_Tr2,Rec", "DP_Tr2,Rec", "RP_Tr2,Rec"
    ) |>
        str_replace_all(c("LP" = "latent_period",
                          "DP" = "detection_period",
                          "RP" = "removal_period"))
    
    # Remove repeated sigma and infrat
    beta_in <- str_subset(parnames, "beta")
    if (beta_in[[1]] == "beta_Tr2") {
        model_pars[c(1, 6)] <- "empty"
    } else {
        model_pars[c(11, 16)] <- "empty"
    }
    
    fes <- expand.grid(sildt1,
                       c("trial", "donor", "txd", "weight", "weight1", "weight2")) |>
        rev() |> apply(1, str_flatten, "_")
    
    plt_names <- c(cov_pars, model_pars, fes)
    
    # Some entries like "trial_s" might be missing
    plt_names[plt_names %notin% parnames] <- "empty"
    
    # This clips any rows or columns that are entirely empty
    plt_mat <- matrix(plt_names, nrow = 5)
    plt_mat <- plt_mat[
        which(apply(plt_mat, 1, any_non_empty)),
        which(apply(plt_mat, 2, any_non_empty))
    ]
    plt_names <- as.character(plt_mat)
    
    if (is_empty(plt_names)) return()
    
    pltlst <- with(plts$plts, mget(plt_names))
    
    plt <- plot_grid(plts$title_plt + theme_classic(),
                     plot_grid(plotlist = pltlst, ncol = nrow(plt_mat)),
                     ncol = 1, rel_heights = c(0.05, 1))
    
    plt_str <- str_glue("{post_dir}/{dataset}-s{scen}-{rep}-posteriors")
    ggsave(str_c(plt_str, ".png"), plt, width = 10, height = 12)
    ggsave(str_c(plt_str, ".pdf"), plt, width = 10, height = 12)
    message(str_glue(" - plotted '{plt_str}'"))
    
    
    # Covariance only
    # Do we want to include the latent_period here?
    plt_names <- cov_pars
    plt_names[plt_names %notin% parnames] <- "empty"
    
    pltlst <- with(plts$plts, mget(plt_names))
    
    plt_mat <- matrix(plt_names, nrow = 5)
    plt_mat <- plt_mat[
        which(apply(plt_mat, 1, any_non_empty)),
        which(apply(plt_mat, 2, any_non_empty))
    ]
    plt_names <- as.character(plt_mat)
    
    if (is_empty(plt_names)) {
        message(" - No Genetic Variance parameters detected")
        return()
    }
    
    cov_plt <- plot_grid(plts$title_plt + theme_classic(),
                         plot_grid(plotlist = pltlst, ncol = 3),
                         ncol = 1, rel_heights = c(0.05, 1))
    
    plt_str <- str_glue("{cov_dir}/{dataset}-s{scen}-{rep}-covs") |>
        str_replace("scen-", "s")
    ggsave(str_c(plt_str, ".png"), cov_plt, width = 9, height = 8)
    ggsave(str_c(plt_str, ".pdf"), cov_plt, width = 8, height = 8)
    message(str_glue(" - Plotted '{plt_str}'"))
    
    plt
}

if (FALSE) {
    # get_posteriors(dataset = "testing", scens = 1:2)
    dataset <- "fb-test"; scens <- 1:4
    walk(scens, \(i) get_posterior(dataset, i, 1))
    
    dataset <- "fb-qtest"; scens <- 1:13
    walk(scens, \(i) get_posterior(dataset, i, 1))
    
    # get_posteriors(dataset = "sim-base-inf", scens = 5, reps = 1:5, combine = FALSE)
    # walk(scens, \(i) get_posterior(dataset, i, 1))
}
