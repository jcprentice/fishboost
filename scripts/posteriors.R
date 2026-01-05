{
    library(data.table)
    library(stringr)
    library(purrr)
    library(cowplot)
    source("plotting/plot_posteriors.R")
}

get_posteriors <- function(dataset = "fb-final", scens = 1:2, reps = 1,
                           combine = FALSE) {
    # dataset <- "fb-test"; scens <- 1:5; reps <- 1:5; combine <- FALSE
    # dataset <- "sim-base-inf"; scens <- 1:10; reps <- 0; combine <- TRUE
    
    base_dir <- str_glue("datasets/{dataset}")
    res_dir  <- str_glue("{base_dir}/results")
    gfx_dir  <- str_glue("{base_dir}/gfx")
    post_dir <- str_glue("{gfx_dir}/posteriors")
    cov_dir  <- str_glue("{gfx_dir}/posteriors_covs")
    
    walk(c(gfx_dir, post_dir, cov_dir), \(d) {
        if (!dir.exists(d)) dir.create(d, recursive = TRUE)
    })
    
    
    rfs <- list.files(res_dir) |>
        str_remove_all("scen-|.rds") |>
        str_sort(numeric = TRUE) |>
        data.table(name = _)
    
    rfs[, `:=`(scen = name |> str_split_i("-", 1) |> as.integer(),
               rep  = name |> str_split_i("-", 2) |> as.integer())]
    
    if (!is.null(scens) && any(scens != 0)) rfs <- rfs[scen %in% scens]
    if (!is.null(reps)  && any(reps != 0))  rfs <- rfs[rep  %in% reps]
    if (combine) rfs <- rfs[, .(name = str_c(scen, "-", 0), rep = 0), scen]
    
    scen_reps <- rfs$name
    
    plts <- map(scen_reps, possibly(~ plot_posteriors(
        dataset = dataset,
        name = .x,
        ci = "hpdi",
        draw = "density")
    ))
    
    
    
    # iwalk(plts, \(x, i) {
    #     if (is.null(x)) return()
    #     
    #     pdf_str <- str_glue("{post_dir}/{dataset}-{i}-posteriors.png")
    #     ggsave(pdf_str, x$pars, width = 10, height = 12)
    #     message(str_glue("plotted '{pdf_str}'"))
    # })
    
    
    # Create a plot with subplots aligned by trait
    empty <- ggplot() + theme_classic()
    sildt1 <- c("s", "i", "l", "d", "t")
    sildt2 <- str_c(sildt1, sildt1)
    any_non_empty <- function(x) any(x != "empty")
    
    walk(seq_along(scen_reps), \(i) {
        x <- plts[[i]]
        if (is.null(x)) return()
        
        x$plts$empty <- empty
        parnames <- names(x$plts)
        
        fes <- expand.grid(sildt1,
                           c("trial", "donor", "txd", "weight", "weight1", "weight2")) |>
            rev() |> apply(1, str_flatten, "_")
        
        # Do we want to include the latent_period here?
        bici_pars <- c(
            "sigma",  "beta_Tr1", "LP_Tr1,Don", "DP_Tr1,Don", "RP_Tr1,Don", 
            "infrat", "empty",    "LP_Tr1,Rec", "DP_Tr1,Rec", "RP_Tr1,Rec",
            "sigma",  "beta_Tr2", "LP_Tr2,Don", "DP_Tr2,Don", "RP_Tr2,Don",
            "infrat", "empty",    "LP_Tr2,Rec", "DP_Tr2,Rec", "RP_Tr2,Rec"
        )
        
        # Remove repeated sigma and infrat
        beta_in <- str_subset(parnames, "beta")
        if (beta_in[[1]] == "beta_Tr2") {
            bici_pars[c(1, 6)] <- "empty"
        } else {
            bici_pars[c(11, 16)] <- "empty"
        }
        
        sire_pars <- c("sigma", "beta", "LP", "DP", "RP")
        
        cov_pars <- c(str_c("cov_G_", sildt2),
                      "r_G_si", "r_G_st", "empty", "empty", "r_G_it",
                      str_c("cov_E_", sildt2),
                      str_c("cov_P_", sildt2))
        plt_names <- c(cov_pars,
                       if (any(str_detect(parnames, "Tr"))) bici_pars else sire_pars,
                       fes) |>
            str_replace_all(c("LP" = "latent_period",
                              "DP" = "detection_period",
                              "RP" = "removal_period"))
        
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
        
        pltlst <- with(x$plts, mget(plt_names))
        
        plt <- plot_grid(x$title_plt + theme_classic(),
                         plot_grid(plotlist = pltlst, ncol = nrow(plt_mat)),
                         ncol = 1, rel_heights = c(0.05, 1))
        
        plt_str <- str_glue("{post_dir}/{dataset}-s{sr}-posteriors",
                            sr = str_remove(scen_reps[[i]], "-0"))
        # ggsave(str_c(plt_str, ".png"), plt, width = 10, height = 12)
        ggsave(str_c(plt_str, ".pdf"), plt, width = 10, height = 12)
        message(str_glue(" - plotted '{plt_str}'"))
    })
    
    walk(seq_along(scen_reps), \(i) {
        x <- plts[[i]]
        if (is.null(x)) return()
        
        x$plts$empty <- empty
        parnames <- names(x$plts)
        
        # Do we want to include the latent_period here?
        plt_names <- c("cov_G_ss", "cov_G_ii", "cov_G_tt",
                       "r_G_si",   "r_G_st",   "r_G_it",
                       "cov_E_ss", "cov_E_ii", "cov_E_tt",
                       "cov_P_ss", "cov_P_ii", "cov_P_tt",
                       "h2_ss",    "h2_ii",    "h2_tt")
        plt_names[plt_names %notin% parnames] <- "empty"
        
        pltlst <- with(x$plts, mget(plt_names))
        
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
        
        plt <- plot_grid(x$title_plt + theme_classic(),
                         plot_grid(plotlist = pltlst, ncol = 3),
                         ncol = 1, rel_heights = c(0.05, 1))
        
        plt_str <- str_glue("{cov_dir}/{dataset}-s{scen_reps[[i]]}-covs") |>
            str_replace("scen-", "s")
        # ggsave(str_c(plt_str, ".png"), plt, width = 9, height = 8)
        ggsave(str_c(plt_str, ".pdf"), plt, width = 8, height = 8)
        message(str_glue(" - Plotted '{plt_str}'"))
    })
}

if (FALSE) {
    # get_posteriors(dataset = "testing", scens = 1:2)
    get_posteriors(dataset = "fb-test", scens = 1:4, reps = 1)
    get_posteriors(dataset = "fb-qtest", scens = 1:13, reps = 1)
    # get_posteriors(dataset = "sim-base-inf", scens = 5, reps = 1:5, combine = FALSE)
}
