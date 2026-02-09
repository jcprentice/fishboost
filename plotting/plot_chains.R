{
    library(data.table)
    library(stringr)
    library(purrr)
    library(ggplot2)
    library(ggeasy)
    library(RColorBrewer)
    library(cowplot)
    
    source("rename_pars.R")
}

plot_chains <- function(dataset = "fb-test", scen = 1, rep = 1) {
    # dataset <- "sim-base-inf"; scen <- 2; rep <- 1
    # dataset <- "fb-test"; scen <- 1; rep <- 1
    message(str_glue("plotting chain for '{dataset}' / s{scen}-{rep}"))
    
    {
        base_dir <- str_glue("datasets/{dataset}")
        data_dir <- str_glue("{base_dir}/data")
        res_dir <- str_glue("{base_dir}/results")
        gfx_dir <- str_glue("{base_dir}/gfx")
        chains_dir <- str_glue("{gfx_dir}/gfx/chains")
    }
    
    c(gfx_dir, chains_dir) |>
        discard(dir.exists) |>
        walk(~ message("- mkdir ", .x)) |>
        walk(dir.create, recursive = TRUE)
    
    trace_files <- list.files(str_glue("{data_dir}/scen-{scen}-{rep}-out/output-inf"),
                              "param_", full.names = TRUE) |>
        str_sort(numeric = TRUE)
    
    traces <- map(trace_files, fread) |>
        rbindlist(idcol = "chain")
    
    traces[, str_subset(names(.SD), "State|Prior|L\\^|^G_|N\\^") := NULL]
    setnames(traces, rename_bici_pars)
    
    ltraces <- traces |> melt(id.vars = "chain",
                              variable.name = "parameter",
                              value.name = "y")
    ltraces[, x := seq(.N), parameter]
    ltraces[, chain := factor(chain)]
    
    pars <- ltraces$parameter |> levels()
    
    rf <- str_glue("{res_dir}/scen-{scen}-{rep}.rds")
    if (file.exists(rf)) {
        res_file <- readRDS(rf)
        priors <- res_file$params$priors
        priors[, parameter := str_replace_all(parameter, c("latent_period" = "LP",
                                                           "detection_period" = "DP",
                                                           "removal_period" = "RP"))]
        priors <- priors[parameter %in% pars,
                         .(parameter, min = val1, max = val2, true_val)]
        nchains <- res_file$params$nchains
        description <- res_file$params$description
    } else {
        priors <- ltraces[, .(min = min(y), max = max(y),
                              true_val = NA_real_),
                          parameter]
        nchains <- ltraces$chain |> levels() |> length()
        description <- ""
    }
    
    thin <- 1L
    plts <- map(pars, \(par) {
        # par <- "cov_G_ii"
        X <- ltraces[parameter == par][seq(1, .N, thin)]
        ytrue <- if (str_detect(dataset, "sim")) {
            priors[parameter == par, true_val]
        }
        
        yrng <- ltraces[parameter == par, range(ytrue, y)]
        yrng <- yrng + c(-0.1, 0.1) * abs(yrng)
        
        ggplot(X, aes(x, y, group = par, colour = chain)) +
            geom_point(size = 0.1) +
            {if (str_detect(dataset, "sim"))
                geom_hline(yintercept = ytrue,
                           linetype = "dashed",
                           colour = "green")} +
            scale_colour_manual(values = rep(brewer.pal(4, "Set1"),
                                             length = nchains)) +
            scale_y_continuous(limits = ~ range(.x, 0, yrng)) +
            labs(x = NULL,
                 y = "value",
                 title = rename_pars(par)) +
            theme_classic() +
            theme(plot.title = element_text(size = 12),
                  plot.subtitle = element_text(size = 10),
                  legend.position = "none",
                  axis.text.x = element_blank())
    }) |>
        setNames(pars)
    
    plts$empty <- ggplot() + theme_classic()
    
    sildt1 <- str_chars("sildt")
    sildt2 <- str_c(sildt1, sildt1)
    any_non_empty <- function(x) any(x != "empty")
    
    cov_pars <- c(str_c("cov_G_", sildt2),
                  "r_G_si", "r_G_st", "empty", "empty", "r_G_it",
                  str_c("cov_E_", sildt2), str_c("cov_P_", sildt2))
    
    model_pars <- c(
        "sigma",  "beta_Tr1", "LP_Tr1,Don", "DP_Tr1,Don", "RP_Tr1,Don",
        "infrat", "empty",    "LP_Tr1,Rec", "DP_Tr1,Rec", "RP_Tr1,Rec",
        "sigma",  "beta_Tr2", "LP_Tr2,Don", "DP_Tr2,Don", "RP_Tr2,Don",
        "infrat", "empty",    "LP_Tr2,Rec", "DP_Tr2,Rec", "RP_Tr2,Rec")
    
    fes <- expand.grid(sildt1,
                       c("trial", "donor", "txd", "weight", "weight1", "weight2")) |>
        rev() |> apply(1, str_flatten, "_")
    
    plt_names <- c(cov_pars, model_pars, fes)
    
    if ("beta_Tr1" %notin% pars) {
        idx <- str_which(plt_names, "sigma")[[1]]
        plt_names[seq(idx, idx + 9)] <- "empty"
    }
    if ("beta_Tr2" %notin% pars) {
        idx <- str_which(plt_names, "sigma")[[2]]
        plt_names[seq(idx, idx + 9)] <- "empty"
    }
    
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
    
    title_plt <- ggplot() +
        labs(title = str_glue("{dataset} / s{scen}-{rep}"),
             subtitle = description) +
        theme_classic()
    theme(plot.title = element_text(size = 12),
          plot.subtitle = element_text(size = 10))
    
    plt <- plot_grid(title_plt,
                     plot_grid(plotlist = pltlst,
                               ncol = nrow(plt_mat),
                               align = "v"),
                     ncol = 1, rel_heights = c(0.08, 1))
    
    if (!dir.exists(gfx_dir)) {
        message("- mkdir ", gfx_dir)
        dir.create(gfx_dir)
    }
    # PDFs are huge here
    ggsave(str_glue("{gfx_dir}/{dataset}-s{scen}-{rep}-chains.png"),
           plt, height = 9, width = 12)
    plt
}

if (FALSE) {
    dataset <- "fb-qtest"
    scens <- str_glue("datasets/{dataset}/results") |>
        list.files() |> str_split_i("-", 2) |> as.integer() |> sort() |> unique()
    
    walk(scens, ~ plot_chains(dataset, .x, 1))
}
