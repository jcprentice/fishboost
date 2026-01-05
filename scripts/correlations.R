{
    library(stringr)
    library(data.table)
    library(purrr)
    library(ggplot2)
    library(ggcorrplot)
}

get_correlations <- function(dataset = "fb-final", scen = 1, rep = 1) {
    # dataset <- "fb-test"; scen <- 1; rep <- 1;

    {
        base_dir <- str_glue("datasets/{dataset}")
        data_dir <- str_glue("{base_dir}/data")
        res_dir  <- str_glue("{base_dir}/results")
        gfx_dir  <- str_glue("{base_dir}/gfx")
    }
    
    x <- if (rep == 0) {
        files <- list.files(data_dir, "trace_combine.tsv",
                   recursive = TRUE, full.names = TRUE) |>
            str_sort(numeric = TRUE) |>
            str_subset(str_glue("scen-{scen}-"))
        if (length(files) == 0) {
            return(NULL)
        }
        files |> map(fread) |> rbindlist()
    } else {
        f <- str_glue("{data_dir}/scen-{scen}-{rep}-out/trace_combine.tsv")
        if (!file.exists(f)) {
            return(NULL)
        }
        fread(f)
    }
    
    params <- if (rep == 0) {
        list.files(res_dir, full.names = TRUE) |>
            str_subset(str_glue("scen-{scen}")) |> str_sort(numeric = TRUE) |>
            first() |> readRDS() |> _$params
    } else {
        str_glue("{res_dir}/scen-{scen}-{rep}.rds") |>
            readRDS() |> _$params
    }

    x[, str_subset(names(x), "Group|state|^G_") := NULL]
    
    setnames(x, str_replace_all(names(x),
                                c("latent_period" = "LP",
                                  "detection_period" = "DP",
                                  "removal_period" = "RP",
                                  "," = "_")))
    nx <- names(x)
    
    # Genetic covariances and correlations
    GVs <- c("cov_G_ss", "cov_G_ii", "cov_G_tt", "r_G_si", "r_G_st", "r_G_it")
    EVs <- str_replace(GVs, "G", "E")
    # parameter associated with each fixed effect
    fe_par <- c(s = "beta", i = "beta", l = "LP", d = "DP", t = "RP")
    # FE names
    fes <- c("trial", "donor", "txd", "weight")
    
    
    sildt <- c("s", "i", "l", "d", "t")
    
    my_theme <- function() {
        theme_classic() +
            theme(axis.title = element_blank(),
                  axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    }
    
    plts <- map(sildt, \(fe) {
        pars <- c(GVs, fe_par[[fe]], str_c(fes, "_", fe))
        
        # Some sets don't have weight
        if (!any(str_detect(nx, "weight"))) {
            pars <- str_subset(pars, "^weight", negate = TRUE)
        }
        
        # Weight may be nested
        if (any(str_detect(nx, "weight1"))) {
            pars <- c(str_subset(pars, "weight", negate = TRUE),
                      str_c("weight", 1:2, "_", fe))
        }
        
        # This just makes sure that pars actually exist
        pars <- intersect(pars, nx)
        if (is_empty(pars)) return(NULL)
        
        ggcorrplot(cor(x[, ..pars]), method = "circle") +
            my_theme()
        
        # datasets/dataset/gfx/dataset-N-cor_l.png
        # ggsave(str_glue("{gfx_dir}/correlations/{dataset}-s{scen}-corr_{fe}.png"),
        #        plt, width = 6, height = 6)
    }) |>
        setNames(sildt)
    
    description <- params$description |>
        str_remove_all(", convergence|, GRM \\w*") |>
        str_replace_all(c("inf_model 1" = "inf: I = D",
                          "inf_model 2" = "inf: I = 0.1*D",
                          "inf_model 3" = "inf: Don = 0.1*Rec",
                          "inf_model 4" = "inf: Don = r*Rec"))
    
    title_str <- str_glue("{dataset}/s{scen}-{rep}: {description}")
    
    plts$all <- ggcorrplot(cor(x),
                           method = "circle",
                           title = title_str) +
        my_theme()
    
    GxEpars <- intersect(c(GVs, EVs), names(x))
    plts$GxE <- if (length(GxEpars) > 0) {
        ggcorrplot(cor(x[, ..GxEpars]),
                   method = "circle",
                   title = title_str) +
            my_theme()
    }
        
    pars_G_lp <- str_subset(names(x), "cov_G|LP|DP|RP")
    plts$G_lp <- if (length(pars_G_lp) > 0) {
        ggcorrplot(cor(x[, ..pars_G_lp]),
                                 method = "circle",
                                 title = title_str) +
        my_theme()
    }
    
    plts
}
# dataset <- "fb-parasites4"; scen <- 3
# dataset <- "fb-fes2"; scens <- 1:11
# dataset <- "fb-fes3"; scens <- 1:11
dataset <- "fb-test"; scens <- 1:4; reps <- 0;

scen_reps <- expand.grid(rep = reps, scen = scens) |> rev() |> as.data.table()
names <- expand.grid(reps, "-", scens, "s") |> rev() |>
    apply(1, str_flatten) |> str_remove_all(" |-0")

plts <- map(seq_len(nrow(scen_reps)), \(i) {
    get_correlations(dataset = dataset,
                     scen = scen_reps[i, scen],
                     rep = scen_reps[i, rep])
}) |>
    setNames(names)

sildt <- c("s", "i", "l", "d", "t")

cor_dir <- str_glue("datasets/{dataset}/gfx/correlations")
if (!dir.exists(cor_dir)) dir.create(cor_dir)

walk(seq_len(nrow(scen_reps)), \(i) {
    scen <- scen_reps[i, scen]
    rep  <- scen_reps[i, rep]
    name <- names[[i]]
    
    
    # walk(sildt, \(fe) {
    #     if (!is.null(plts[[name]][[fe]])) {
    #         ggsave(str_glue("{cor_dir}/{dataset}-{name}-corr_{fe}.png"),
    #                plts[[name]][[fe]], width = 6, height = 6)
    #     }
    # })
    
    ggsave(str_glue("{cor_dir}/{dataset}-{name}-corr_all.pdf"),
           plts[[name]]$all, width = 12, height = 12)
    
    if (!is.null(plts[[name]]$GxE)) {
        ggsave(str_glue("{cor_dir}/{dataset}-{name}-corr_GxE.pdf"),
               plts[[name]]$GxE, width = 6, height = 6)
    }
    
    if (!is.null(plts[[name]]$G_lp)) {
        ggsave(str_glue("{cor_dir}/{dataset}-{name}-corr_G_lp.pdf"),
               plts[[name]]$G_lp, width = 12, height = 12)
    }
})
