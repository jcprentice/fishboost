{
    library(data.table)
    library(purrr)
    library(stringr)
    library(ggplot2)
    library(cowplot)
    
    source("rename_pars.R")
}

pars_boxplots <- function(dataset = "fb-final", scens = 0, st_str = "") {
    # dataset <- "fb-lp3"; scens <- 1:12
    # st_str = ""
    
    data_dir <- str_glue("datasets/{dataset}/data")
    res_dir  <- str_glue("datasets/{dataset}/results")
    gfx_dir  <- str_glue("datasets/{dataset}/gfx")
    
    files <- list.files(data_dir, "^trace_combine",
                        recursive = TRUE, full.names = TRUE) |>
        str_sort(numeric = TRUE)
    
    if (all(scens != 0)) {
        files <- files |> 
            keep(~ .x |> str_split_i("/", 4) |> str_split_i("-", 2) |>
                     as.integer() |> is.element(scens))
    }
    
    x <- map(files, ~ {
        tc <- fread(.x)
        pars <- names(tc) |>
            str_subset("state|State|^G_|Group effect", negate = TRUE)
        tc[, setdiff(names(tc), pars) := NULL]
        tc[, map(.SD, as.numeric)]
    }) |>
        rbindlist(idcol = "scen", fill = TRUE) |>
        melt(id.vars = "scen",
             variable.name = "parameter",
             variable.factor = FALSE) |>
        split(by = "parameter") |>
        map(~ {
            .x[, `:=`(scen = ordered(str_c("s", scens[scen]),
                                     levels = rev(str_c("s", unique(scens)))),
                      parameter = NULL)]
        })
    
    pars <- names(x)
    tidy_pars <- setNames(rename_pars(pars), pars)
    
    f <- list.files(str_glue("{res_dir}"), "scen", full.names = TRUE)[[1]]
    priors <- readRDS(f)$params$priors
    
    lims <- map(x, ~ .x[, .(min_val = min(value), max_val = max(value))]) |>
        rbindlist(idcol = "parameter")
    
    plts <- map(pars, \(par) {
        # par <- "beta"
        xp <- x[[par]]
        rng <- lims[parameter == par,
                    .(min = min(c(min_val, xp$value)),
                      max = max(c(min_val, xp$value)))] |> as.list()
        
        ggplot(xp, aes(x = scen, y = value)) +
            geom_boxplot(fill = "tomato",
                         outliers = FALSE,
                         staplewidth = 1) +
            expand_limits(y = 0) +
            ylim(rng$min, rng$max) +
            # coord_flip() +
            labs(x = "Scenario",
                 y = "Value",
                 title = tidy_pars[[par]]) +
            theme_classic()
    }) |> setNames(pars)
    
    title_plt <- ggplot() +
        labs(title = str_glue("Dataset: '{dataset}'"),
             subtitle = st_str) +
        theme_classic() +
        theme(plot.title = element_text(size = 22),
              plot.subtitle = element_text(size = 16))
    
    plts$empty <- ggplot() + theme_classic()
    
    sildt1 <- c("s", "i", "l", "d", "t")
    sildt2 <- str_c(sildt1, sildt1)
    any_non_empty <- function(x) any(x != "empty")
    
    cov_pars <- c(str_c("cov_G_", sildt2),
                  "r_G_si", "r_G_st", "empty", "empty", "r_G_it",
                  str_c("cov_E_", sildt2),
                  str_c("cov_P_", sildt2))
    
    model_pars <- c(
        "sigma", "beta_Tr1", "LP_Tr1,Don", "DP_Tr1,Don", "RP_Tr1,Don", 
        "empty", "empty",    "LP_Tr1,Rec", "DP_Tr1,Rec", "RP_Tr1,Rec",
        "empty", "beta_Tr2", "LP_Tr2,Don", "DP_Tr2,Don", "RP_Tr2,Don",
        "empty", "empty",    "LP_Tr2,Rec", "DP_Tr2,Rec", "RP_Tr2,Rec"
    ) |>
        str_replace_all(c("LP" = "latent_period",
                          "DP" = "detection_period",
                          "RP" = "removal_period"))
    
    fes <- expand.grid(sildt1,
                       c("trial", "donor", "txd", "weight", "weight1", "weight2")) |>
        rev() |> apply(1, str_flatten, "_")
    
    plt_names <- c(cov_pars, model_pars, fes)
    
    # Some entries like "trial_s" might be missing
    plt_names[plt_names %notin% pars] <- "empty"
    
    # This clips any rows or columns that are entirely empty
    plt_mat <- matrix(plt_names, nrow = 5)
    plt_mat <- plt_mat[
        which(apply(plt_mat, 1, any_non_empty)),
        which(apply(plt_mat, 2, any_non_empty))
    ]
    plt_names <- c(plt_mat)
    
    pltlst <- with(plts, mget(plt_names))
    
    plt <- plot_grid(title_plt,
                     plot_grid(plotlist = pltlst, ncol = nrow(plt_mat)),
                     ncol = 1, rel_heights = c(0.06, 1))
    
    ggsave(str_glue("{gfx_dir}/{dataset}-all_boxplots.png"), plt,
           width = 15, height = 18)
    ggsave(str_glue("{gfx_dir}/{dataset}-all_boxplots.pdf"), plt,
           width = 15, height = 18)
}

if (FALSE) {
    pars_boxplots("fb-test", 1:5, "Testing BICI on FB data")
    # pars_boxplots("fb-final", 1:8, "Testing Unlinked vs Linked Traits and FEs")
    # pars_boxplots("fb-final2", 1:4, "Testing Unlinked vs Linked Traits and FEs")
    # pars_boxplots("fb-lp", 1:12, "Testing varying the LP")
    # pars_boxplots("fb-donors", 1:3, "Testing reclassifying Seeder fish")
    # pars_boxplots("fb-simple", 1:6, "Testing if we can fit the LP")
    # pars_boxplots("fb-simple-b", 1:6, "Testing if we can fit the LP with BICI")
}
