library(pROC)
library(ggplot2)
library(cowplot)

mark_cutoff <- function(ranks, p = 0.8) {
    ranks2 <- copy(ranks)
    cutoff <- p * ranks2[, length(unique(sus))]
    ranks2[, `:=`(
        psus = factor(fifelse(sus > cutoff, 1L, 0L), c(0L, 1L)),
        pinf = factor(fifelse(inf > cutoff, 1L, 0L), c(0L, 1L)),
        prec = factor(fifelse(rec > cutoff, 1L, 0L), c(0L, 1L))
    )]

    ranks2
}

plot_roc_curves <- function(ranks, p = 0.8) {

    ranks2 <- mark_cutoff(ranks, p)

    roc_sus <- roc(ranks2$psus, ranks2$est_sus)
    roc_inf <- roc(ranks2$pinf, ranks2$est_inf)
    roc_rec <- roc(ranks2$prec, ranks2$est_rec)

    auc_sus <- signif(roc_sus$auc, digits = 3L)
    auc_inf <- signif(roc_inf$auc, digits = 3L)
    auc_rec <- signif(roc_rec$auc, digits = 3L)

    cat(paste0("AUC for susceptibility = ", auc_sus, "\n"))
    cat(paste0("AUC for infectivity = ",    auc_inf, "\n"))
    cat(paste0("AUC for recoverability = ", auc_rec, "\n"))

    plt_roc <- ggroc(list(roc_sus, roc_inf, roc_rec)) +
        geom_line(size = 1.2) +
        labs(title = "ROC",
             x = "Specificity",
             y = "Sensitivity") +
        scale_colour_manual(labels = c(paste0("Susceptibility (AUC = ", auc_sus, ")"),
                                       paste0("Infectivity (AUC = ", auc_inf, ")"),
                                       paste0("Recoverability (AUC = ", auc_rec, ")")),
                            values = c("red", "darkgreen", "blue")) +
        theme(legend.position = c(0.8, 0.3))

    # plt_sus <- ggroc(roc_sus) +
    #     labs(title = paste0("AUC for susceptibility = ", auc_sus),
    #          x = "Specificity",
    #          y = "Sensitivity")
    #
    # plt_inf <- ggroc(roc_inf) +
    #     labs(title = paste0("AUC for infectivity = ", auc_inf),
    #          x = "Specificity",
    #          y = "Sensitivity")
    #
    # plt_rec <- ggroc(roc_rec) +
    #     labs(title = paste0("AUC for recoverability = ", auc_rec),
    #          x = "Specificity",
    #          y = "Sensitivity")
    #
    # plt_roc <- plot_grid(
    #     plotlist = list(plt_sus, plt_inf, plt_rec)
    # )

    plt_roc
}
