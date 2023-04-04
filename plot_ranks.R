{
    library(data.table)
    library(glue)
    library(gtools)
    library(lubridate)
    library(ggplot2)
    library(cowplot)
    library(extrafont)
    
    source("get_ranks.R")
    source("get_roc.R")
    source("utils.R")
}

plot_ranks <- function(data_set = "sim-Gcov", scenario = 1) {
    
    # data_set = "sim-Gsir_cov0_Dlidr-1-mpi"; scenario = 1
    
    if (startsWith(data_set, "fb")) {
        message("Can't get ranks for FB data!")
        return()
    }
    
    # Directory to save all files to
    gfx_dir <- glue("gfx/{data_set}")
    if (!dir.exists(gfx_dir)) {
        dir.create(gfx_dir)
    }
    
    message(glue("Generating rank plots for {data_set}, scenario {scenario} ..."))
    
    BVs <- data.table()
    Rks <- data.table()
    
    # get list of files to scan
    res_dir <- glue("results/{data_set}")
    res_files <- list.files(res_dir, pattern = glue("scen-{scenario}-"))
    
    if (length(res_files) == 0L) {
        message(glue("No files available for Scenario {scenario}!"))
        return(NULL)
    } else {
        message(glue("Found {length(res_files)} files"))
    }
    
    for (i in seq_along(res_files)) {
        # Loads: params, pop, parameter_estimates, estimated_BVs, pred_accs,
        # ranks, time_taken
        f <- glue("{res_dir}/{res_files[i]}")
        load(f)
        
        if (!is.null(estimated_BVs)) {
            ranks <- get_ranks(pop, estimated_BVs, params)
        } else {
            message(glue("{data_set} has no traits"))
            return(NULL)
        }
        
        Rks <- rbind(Rks, ranks)
        
        sires <- pop[sdp == "sire", .I]
        
        BV = data.table()
        if (grepl("s", params$use_traits)) {
            BV <- rbind(BV, cbind(pop[sires, .(id, variable = "sus", true = rank(susceptibility_BV))],
                                  estimated_BVs[sires, .(estimated = rank(susceptibility))]))
        }
        if (grepl("i", params$use_traits)) {
            BV <- rbind(BV, cbind(pop[sires, .(id, variable = "inf", true = rank(infectivity_BV))],
                                  estimated_BVs[sires, .(estimated = rank(infectivity))]))
        }
        if (grepl("r", params$use_traits)) {
            BV <- rbind(BV, cbind(pop[sires, .(id, variable = "rec", true = rank(recoverability_BV))],
                                  estimated_BVs[sires, .(estimated = rank(recoverability))]))
        }
        
        BVs <- rbind(BVs, BV)
    }
    
    # Will reuse this several times
    if (is.null(params$label)) {
        params$label <- scenario
    }
    
    title_plt <- ggplot() + labs(title = glue("Scenario {params$label}, {params$description}"))
    
    
    # Check ranks and ROC curves ----
    if (length(unique(Rks$sus)) <= 1) {
        message("no ranks")
        plt_rank <- NA
    } else {
        message("ranks")
        # Assign true/false positive/negative to all points
        point_colour <- function(true, estimated, cutoff) {
            fifelse(estimated > cutoff & true > cutoff, "TP",
                    fifelse(estimated <= cutoff & true <= cutoff, "TN",
                            fifelse(estimated > cutoff & true <= cutoff, "FP", "FN")))
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
        
        p_rec <- ggplot(BVs[variable == "rec"], aes(x = true, y = estimated)) +
            geom_jitter(aes(x = true, y = estimated, colour = Top20), size = 2.0) +
            geom_smooth(aes(x = true, y = estimated), method = "lm") +
            labs(title = "Recoverability", x = "True rank", y = "Estimated rank")
        
        plt_rank <- plot_grid(
            title_plt,
            plot_grid(plotlist = list(p_sus, p_inf, p_rec, plt_roc)),
            ncol = 1, rel_heights = c(0.06, 1))
        
        ranks_str <- glue("{gfx_dir}/ranks-scen{scenario}.png")
        message("Plotting Ranks: ", ranks_str)
        ggsave(ranks_str, plt_rank, width = 1000, height = 1000, units = "px")
    }
    
    plt_rank
}

