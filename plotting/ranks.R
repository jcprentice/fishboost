{
    library(data.table)
    library(stringr)
    library(gtools)
    library(lubridate)
    library(ggplot2)
    library(cowplot)
    library(extrafont)
    
    source("get_ranks.R")
    source("get_roc.R")
    source("utils.R")
}

plot_ranks <- function(dataset = "sim-Gcov", scen = 1) {
    
    # dataset = "sim-Gsit_cov0_Dildt-1-mpi"; scen = 1
    
    if (str_starts(dataset, "fb")) {
        message("Can't get ranks for FB data!")
        return()
    }
    
    # Directory to save all files to
    gfx_dir <- str_glue("datasets/{dataset}/gfx")
    if (!dir.exists(gfx_dir)) {
        message("- mkdir ", gfx_dir)
        dir.create(gfx_dir)
    }
    
    message(str_glue("Generating rank plots for {dataset}, scenario {scen} ..."))
    
    BVs <- data.table()
    Rks <- data.table()
    
    # get list of files to scan
    res_dir <- str_glue("datasets/{dataset}/results")
    res_files <- list.files(res_dir, pattern = str_glue("scen-{scen}-"))
    
    if (is_empty(res_files)) {
        message(str_glue("No files available for Scenario {scen}!"))
        return(NULL)
    } else {
        message(str_glue("Found {length(res_files)} files"))
    }
    
    for (i in seq_along(res_files)) {
        # Loads: params, popn, parameter_estimates, estimated_BVs, pred_accs,
        # ranks, time_taken
        f <- str_glue("{res_dir}/{res_files[[i]]}")
        load(f)
        
        if (!is.null(estimated_BVs)) {
            ranks <- get_ranks(popn, estimated_BVs, params)
        } else {
            message(str_glue("{dataset} has no traits"))
            return(NULL)
        }
        
        Rks <- rbind(Rks, ranks)
        
        sires <- popn[sdp == "sire", .I]
        
        BV = data.table()
        if (str_detect(params$use_traits, "s")) {
            BV <- rbind(BV, cbind(popn[sires, .(id, variable = "sus", true = rank(sus_g))],
                                  estimated_BVs[sires, .(estimated = rank(sus))]))
        }
        if (str_detect(params$use_traits, "i")) {
            BV <- rbind(BV, cbind(popn[sires, .(id, variable = "inf", true = rank(inf_g))],
                                  estimated_BVs[sires, .(estimated = rank(inf))]))
        }
        if (str_detect(params$use_traits, "t")) {
            BV <- rbind(BV, cbind(popn[sires, .(id, variable = "tol", true = rank(tol_g))],
                                  estimated_BVs[sires, .(estimated = rank(tol))]))
        }
        
        BVs <- rbind(BVs, BV)
    }
    
    # Will reuse this several times
    params$label <- params$label %||% scen
    
    title_plt <- ggplot() + labs(title = str_glue("Scenario {params$label}, {params$description}"))
    
    
    # Check ranks and ROC curves ----
    if (length(unique(Rks$sus)) <= 1) {
        message("no ranks")
        plt_rank <- NA
    } else {
        message("ranks")
        # Assign true/false positive/negative to all points
        point_colour <- function(true, estimated, cutoff) {
            fcase(estimated > cutoff & true > cutoff,   "TP",
                  estimated <= cutoff & true <= cutoff, "TN",
                  estimated > cutoff & true <= cutoff,  "FP",
                  default =                             "FN")
        }
        
        # generate ROC curve
        plt_roc <- NULL # plot_roc_curves(Rks, p = 0.8)
        
        # label correctly ranked
        cutoff <- 0.8 * max(BVs$true)
        BVs[, Top20 := factor(point_colour(true, estimated, cutoff), levels = c("TP", "TN", "FP", "FN"))]
        
        p_sus <- ggplot(BVs[variable == "sus"]) +
            geom_jitter(aes(x = true, y = estimated, colour = Top20), size = 2.0) +
            geom_smooth(aes(x = true, y = estimated), method = "lm") +
            labs(title = "Susceptibility", x = "True rank", y = "Estimated rank")
        
        p_inf <- ggplot(BVs[variable == "inf"], aes(x = true, y = estimated)) +
            geom_jitter(aes(x = true, y = estimated, colour = Top20), size = 2.0) +
            geom_smooth(aes(x = true, y = estimated), method = "lm") +
            labs(title = "Infectivity", x = "True rank", y = "Estimated rank")
        
        p_tol <- ggplot(BVs[variable == "tol"], aes(x = true, y = estimated)) +
            geom_jitter(aes(x = true, y = estimated, colour = Top20), size = 2.0) +
            geom_smooth(aes(x = true, y = estimated), method = "lm") +
            labs(title = "Tolerance", x = "True rank", y = "Estimated rank")
        
        plt_rank <- plot_grid(
            title_plt,
            plot_grid(plotlist = list(p_sus, p_inf, p_tol, plt_roc)),
            ncol = 1, rel_heights = c(0.06, 1))
        
        ranks_str <- str_glue("{gfx_dir}/ranks-s{scen}.png")
        message("Plotting Ranks: ", ranks_str)
        ggsave(ranks_str, plt_rank, width = 6, height = 6)
    }
    
    plt_rank
}

