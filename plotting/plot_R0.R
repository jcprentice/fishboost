{
    library(data.table)
    library(purrr)
    library(stringr)
    library(ggplot2)
    library(scales)
}

plot_R0 <- function(dataset = "fb-final",
                    scen = 1,
                    plot_boxplot = FALSE) {
    
    # dataset <- "fb-final"; scen <- 1; plot_boxplot = TRUE
    
    x <- fread(str_glue("{dataset}/data/scen-{scen}-1-out/trace_combine.tsv"))
    
    x[, str_subset(names(x), "beta|period|_[GE]_", negate = TRUE) := NULL]
    cols <- names(x)
    
    Vcols <- str_subset(cols, "G|E")
    Vcols2 <- str_replace(Vcols, "cov_|r_", "V")
    setnames(x, Vcols, Vcols2)
    
    # Convert correlation to covariance
    x[, `:=`(
        VG_si = VG_si * sqrt(VG_ss * VG_ii),
        VG_st = VG_st * sqrt(VG_ss * VG_tt),
        VG_it = VG_it * sqrt(VG_ii * VG_tt),
        VE_si = VE_si * sqrt(VE_ss * VE_ii),
        VE_st = VE_st * sqrt(VE_ss * VE_tt),
        VE_it = VE_it * sqrt(VE_ii * VE_tt),
        V_ss = VG_ss + VE_ss,
        V_ii = VG_ii + VE_ii,
        V_tt = VG_tt + VE_tt,
        V_si = VG_si + VE_si,
        V_st = VG_st + VE_st,
        V_it = VG_it + VE_it
    )]
    
    x[, r0_base := beta * (detection_period + removal_period)]
    x[, `:=`(
        correction_GT = exp((VG_ss + VG_ii + VG_tt) / 2 + VG_si - VG_st - VG_it),
        # correction_GT = exp(VG_si + VG_st + VG_it),
        # correction_S  = exp(VG_ss / 2),
        # correction_I  = exp(VG_ii / 2),
        # correction_R  = exp(VG_tt / 2),
        # correction_PT = exp((V_ss + V_ii + V_tt) / 2 + V_si - V_st - V_it)
        correction_PT = exp(V_si + V_st + V_it)
    )]
    
    # am I doing this correctly?
    x[, r0_base_sel := beta * exp(-sqrt(mean(VG_ss)) - sqrt(mean(VG_ii))) *
          (detection_period + removal_period * exp(-sqrt(mean(VG_tt))))]
    
    x[, `:=`(
        r0_full_GT = r0_base * correction_GT,
        # r0_full_S  = r0_base * correction_S,
        # r0_full_I  = r0_base * correction_I,
        # r0_full_R  = r0_base * correction_R,
        r0_full_PT = r0_base * correction_PT,
        
        r0_full_GT_sel = r0_base_sel * correction_GT,
        # r0_full_S_sel  = r0_base_sel * correction_S,
        # r0_full_I_sel  = r0_base_sel * correction_I,
        # r0_full_R_sel  = r0_base_sel * correction_R,
        r0_full_PT_sel = r0_base_sel * correction_PT
    )]
    
    xm <- map(x, mean)
    ym <- with(xm, {
        r0_base    <- beta * (detection_period + removal_period)
        r0_full_GT <- r0_base * exp((VG_ss + VG_ii + VG_tt) / 2 + VG_si - VG_st - VG_it)
        r0_full_PT <- r0_base * exp(V_si + V_st + V_it)
        data.table(
            param = c("r0_base", "r0_full_GT", "r0_full_PT"),
            x = c(r0_base, r0_full_GT, r0_full_PT)
        )
    })
    
    r0cols <- str_subset(names(x), "r0")
    y <- melt(x[, ..r0cols],
              measure.vars = r0cols,
              variable.factor = FALSE)
    
    y[, selection := ifelse(str_ends(variable, "_sel"),
                            "one sigma", "none")]
    y[, variable := as.factor(str_remove(variable, "_sel"))]
    
    
    print(y[, .(mean_R0 = mean(value)), .(variable, selection)])
    
    # dys <- split(y, by = c("variable", "selection")) |>
    #     map(\(x) density(x$value))
    # 
    # ym <- data.table(
    #     variable = names(dys),
    #     mean = y[, mean(value), variable]$V1,
    #     median = y[, median(value), variable]$V1,
    #     x = sapply(dys, \(x) x$x[which.max(x$y)]),
    #     y = sapply(dys, \(x) x$y[which.max(x$y)])
    # )
    
    ylabs <- c(bquote(beta / gamma),
               "Bijma (genotype)",
               # "Bijma (S genotype only)",
               # "Bijma (I only)",
               # "Bijma (R only)",
               "Chris (phenotype)")
    
    ggplot(y[selection == "none"],
                  aes(x = value,
                  fill = variable)) +
        # geom_histogram(aes(x = value, fill = variable),
        #                position = "identity",
        #                alpha = 0.5,
        #                breaks = 10^(seq(0.5, 3, 0.01))) +
        geom_density(alpha = 0.5,
                     # bounds = c(-0.5, Inf),
                     trim = FALSE) +
        geom_vline(data = ym,
                   aes(xintercept = x, colour = param),
                   show.legend = FALSE) +
        geom_vline(xintercept = 1,
                   linetype = "dashed") +
        scale_fill_discrete(bquote(R[0]),
                            labels = ylabs) +
        guides(fill = guide_legend(position = "inside")) +
        # geom_label(data = ym,
        #            aes(x = x, y = y, label = signif(x, 3))) +
        scale_x_continuous(trans = "log10",
                           breaks = trans_breaks("log10", \(x) 10^x),
                           labels = trans_format("log10", math_format(10^.x))) +
        # xlim(0, 10) +
        labs(x = bquote(R[0]),
             y = "Density",
             title = bquote("Posterior distribution for" ~ R[0])) +
        # facet_wrap(. ~ selection, nrow = 2,
        #            labeller = label_both) +
        theme_classic() +
        theme(legend.position.inside = c(0.8, 0.8))
    
    filename <- str_glue("{dataset}/gfx/R0/{dataset}-s{scen}-density-R0.pdf")
    ggsave(filename, plot = last_plot(), width = 8, height = 6)
    
    
    if (plot_boxplot) {
        plt2 <- ggplot(y) +
            geom_boxplot(aes(x = value, y = variable, fill = variable),
                         staplewidth = 0.5,
                         outliers = FALSE) +
            scale_y_discrete(labels = ylabs) +
            labs(x = bquote(R[0]),
                 y = NULL,
                 title = bquote("Posterior values for" ~ R[0])) +
            scale_x_continuous(trans = "log10",
                               breaks = trans_breaks("log10", \(x) 10^x),
                               labels = trans_format("log10", math_format(10^.x))) +
            facet_wrap(. ~ selection, nrow = 2,
                       labeller = label_both) +
            theme(legend.position = "none")
        
        ggsave(str_glue("{dataset}/gfx/R0/{dataset}-s{scen}-boxplot-R0.pdf"),
               plot = plt2, width = 8, height = 6)
    }
    
    plt
}

plt1 = plot_R0(dataset = "fb-final",
               scen = 1,
               plot_boxplot = FALSE)
print(plt1)
# plt2 = plot_R0(dataset = "fb-fes3", scen = 4)
# plt3 = plot_R0(dataset = "fb-fes3", scen = 12)


