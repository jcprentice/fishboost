library(pROC)
library(ggplot2)
library(cowplot)

mark_cutoff <- function(ranks, p = 0.8) {
    ranks2 <- copy(ranks)
    cutoff <- p * ranks2[, length(unique(sus))]
    ranks2[, `:=`(
        psus = factor(fifelse(sus > cutoff, 1L, 0L), c(0L, 1L)),
        pinf = factor(fifelse(inf > cutoff, 1L, 0L), c(0L, 1L)),
        ptol = factor(fifelse(tol > cutoff, 1L, 0L), c(0L, 1L))
    )]

    ranks2
}

plot_roc_curves <- function(ranks, p = 0.8) {

    ranks2 <- mark_cutoff(ranks, p)

    roc_sus <- roc(ranks2$psus, ranks2$est_sus)
    roc_inf <- roc(ranks2$pinf, ranks2$est_inf)
    roc_tol <- roc(ranks2$ptol, ranks2$est_tol)

    auc_sus <- signif(roc_sus$auc, digits = 3L)
    auc_inf <- signif(roc_inf$auc, digits = 3L)
    auc_tol <- signif(roc_tol$auc, digits = 3L)

    cat(str_c("AUC for susceptibility = ", auc_sus, "\n"))
    cat(str_c("AUC for infectivity = ",    auc_inf, "\n"))
    cat(str_c("AUC for tolerance = ",      auc_tol, "\n"))

    plt_roc <- ggroc(list(roc_sus, roc_inf, roc_tol)) +
        geom_line(size = 1.2) +
        labs(title = "ROC",
             x = "Specificity",
             y = "Sensitivity") +
        scale_colour_manual(labels = c(str_c("Susceptibility (AUC = ", auc_sus, ")"),
                                       str_c("Infectivity (AUC = ", auc_inf, ")"),
                                       str_c("Tolerance (AUC = ", auc_tol, ")")),
                            values = c("red", "darkgreen", "blue")) +
        theme(legend.position = c(0.8, 0.3))

    # plt_sus <- ggroc(roc_sus) +
    #     labs(title = str_c("AUC for susceptibility = ", auc_sus),
    #          x = "Specificity",
    #          y = "Sensitivity")
    #
    # plt_inf <- ggroc(roc_inf) +
    #     labs(title = str_c("AUC for infectivity = ", auc_inf),
    #          x = "Specificity",
    #          y = "Sensitivity")
    #
    # plt_tol <- ggroc(roc_tol) +
    #     labs(title = str_c("AUC for tolerance = ", auc_tol),
    #          x = "Specificity",
    #          y = "Sensitivity")
    #
    # plt_roc <- plot_grid(
    #     plotlist = list(plt_sus, plt_inf, plt_tol)
    # )

    plt_roc
}
