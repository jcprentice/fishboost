{
    library(data.table)
    library(stringr)
    library(purrr)
    library(HDInterval)
    library(ggplot2)
    library(ggtext)
    library(ggeasy)
    library(cowplot)
    
    source("rename_pars.R")
}

# This are for calculating top x% over bottom x%
sigma2f <- function(sigma = 0.5, mu = 0, q = 0.9) {
    z <- qnorm(q)
    t <- exp(-z^2 / 2) / (sqrt(2 * pi) * (1 - q))
    topq <- sigma * t + mu
    exp(2 * topq)
}


plot_covariances <- function(dataset = "fb-final",
                             scen = 1L,
                             itn = 1L,
                             ci = "hpdi") {
    
    # dataset <- "fb-parasites"; scen <- 1; itn <- 1; ci <- "hpdi"
    # dataset <- "fb-final"; scen <- 1; itn <- 1; ci <- "hpdi"
    
    base_dir <- str_glue("datasets/{dataset}")
    data_dir <- str_glue("{base_dir}/data")
    res_dir  <- str_glue("{base_dir}/results")
    gfx_dir  <- str_glue("{base_dir}/gfx")
    
    title_plot_str <- str_glue("{dataset} / s{scen}")
    
    res_file <- str_glue("{res_dir}/scen-{scen}-{itn}.rds")
    if (file.exists(res_file)) {
        params <- readRDS(res_file)$params
        title_plot_str <- str_glue("{title_plot_str}: {params$description}")
        burn_prop <- params$burnprop
    } else {
        message(str_glue("- results file '{res_file}' not found, continuing without it"))
        burn_prop <- 0.2
    }
    
    trace_file <- str_glue("{data_dir}/scen-{scen}-{itn}-out/trace_combine.tsv")
    if (!file.exists(trace_file)) {
        message("- No trace file, exiting.")
        return(NULL)
    }
    x <- fread(trace_file)
    x[, str_subset(names(x), "_[GE]_", negate = TRUE) := NULL]
    
    x[, `:=`(cov_P_ss = cov_G_ss + cov_E_ss,
             cov_P_ii = cov_G_ii + cov_E_ii,
             cov_P_tt = cov_G_tt + cov_E_tt)]
    x[, `:=`(h2_ss = cov_G_ss / cov_P_ss,
             h2_ii = cov_G_ii / cov_P_ii,
             h2_tt = cov_G_tt / cov_P_tt)]
    
    # Names as in trace
    pars <- names(x)
    
    # Names as wanted on figure
    pars2 <- rename_pars(pars) |>
        str_replace_all(c("h\\^2 " = "h<sup>2</sup>",
                          " ([G|E|P]) " = "<sub>\\1</sub>")) |>
        setNames(pars) |>
        as.list()
    
    # x2 is x in tidy format
    x2 <- x |>
        melt(measure.vars = pars,
             variable.name = "parameter",
             value.name = "val")
    
    # x3 is parameter / median (currently unused)
    x3 <- x2[, .(median = median(val)), parameter]
    
    
    plts <- map(pars, \(par) {
        # par <- pars[[1]]
        param2 <- pars2[[par]]
        
        xp_val <- x2[parameter == par, val]
        
        dens <- density(xp_val, adjust = 0.5, cut = 0)
        dd <- with(dens, data.table(x, y))
        
        sig_figs <- 3
        mx <- median(xp_val)
        if (ci == "hpdi") {
            hpdi <- hdi(xp_val, credMass = 0.95)
            lq <- hpdi[["lower"]]
            uq <- hpdi[["upper"]]
            # list(lq = lq, uq = uq)
        } else {
            lq <- quantile(xp_val, 0.025)
            uq <- quantile(xp_val, 0.975)
        }
        
        if (par == "cov_G_ii") {
            # rng <- qnorm(c(0.1, 0.9), 0, sqrt(mx)) |> diff() |> exp()
            rng <- sigma2f(sigma = mx, mu = 0, q = 0.9)
            rng_str <- str_c(" x", signif(rng, 3))
            print(str_glue("G cov(inf) = {signif(mx, 3)} ({signif(lq, 3)}, {signif(uq, 3)}), {rng_str}"))
        }
        
        h2_str <- ""
        if (str_starts(par, "h2")) {
            message(str_glue("{par} in ({signif(lq, 3)}, {signif(uq, 3)})"))
            h2 <- x[, mean(get(par))]
            # h2_str <- str_glue(", h<sup>2</sup>={signif(h2, 2)}")
        }
        
        
        if (str_starts(par, "cov")) {
            x_min <- 0; x_max <- 1.1 * uq
        } else if (str_starts(par, "h2")) {
            x_min <- 0; x_max <- 1
        } else if (str_starts(par, "r_")) {
            x_min <- -1; x_max <- 1
        } else {
            stop(str_glue("Bad parameter '{par}'"))
        }
        
        
        ggplot() +
            geom_line(data = dd, aes(x = x, y = y)) +
            geom_area(data = subset(dd, lq < x & x < uq),
                      aes(x = x, y = y),
                      fill = "red", colour = NA, alpha = 0.5) +
            geom_vline(xintercept = mx, linewidth = 1, colour = "blue") +
            coord_cartesian(xlim = c(x_min, x_max)) +
            labs(x = "Value",
                 y = "Density",
                 title = str_glue("{param2}: {round(mx, 2)} ({round(lq, 2)}, {round(uq, 2)}){h2_str}")) +
            theme_classic() +
            theme(plot.title = element_markdown(),
                  text = element_text(size = 10),
                  axis.text.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  # text = element_text(size = 28),
                  legend.position = "none")
    }) |>
        setNames(pars)
    
    plts$empty = ggplot() + theme_classic()
    
    pltG_list <- c("cov_G_ss", "r_G_si",   "r_G_st",
                   "empty",    "cov_G_ii", "r_G_it",
                   "empty",    "empty",    "cov_G_tt")
    pltE_list <- str_replace(pltG_list, "G", "E")
    pltH_list <- str_replace(pltG_list, "cov_G_", "h2_")
    pltP_list <- c("cov_P_ss", "empty",    "empty",
                   "empty",    "cov_P_ii", "empty",
                   "empty",    "empty",    "cov_P_tt")
    
    pltG <- plot_grid(plotlist = plts[pltG_list])
    pltE <- plot_grid(plotlist = plts[pltE_list])
    pltP <- plot_grid(plotlist = plts[pltP_list])
    pltH <- plot_grid(plotlist = plts[pltH_list])
    
    title_plt <- ggplot() + labs(title = title_plot_str) + theme_classic()
    
    pltG_cov <- plot_grid(title_plt, pltG, ncol = 1, rel_heights = c(0.08, 1))
    pltE_cov <- plot_grid(title_plt, pltE, ncol = 1, rel_heights = c(0.08, 1))
    pltP_cov <- plot_grid(title_plt, pltP, ncol = 1, rel_heights = c(0.08, 1))
    pltH_cov <- plot_grid(title_plt, pltH, ncol = 1, rel_heights = c(0.08, 1))
    pltGE_cov <- plot_grid(title_plt, plot_grid(pltG, pltE, ncol = 2), ncol = 1, rel_heights = c(0.08, 1))
    # pltG_cov <- pltG
    # pltE_cov <- pltE
    # pltP_cov <- pltP
    # pltH_cov <- pltH
    # pltGE_cov <- plot_grid(pltG, pltE, ncol = 2)
    
    
    cov_dir <- str_glue("{gfx_dir}/cov")
    if (!dir.exists(cov_dir)) {
        message("- mkdir ", cov_dir)
        dir.create(cov_dir, recursive = TRUE)
    }
    
    walk(c("G", "E", "P", "H", "GE"), \(x) {
        png_str <- str_glue("{cov_dir}/{dataset}-s{scen}-{itn}-cov-{x}.png")
        message(str_glue("Plotting '{png_str}'"))
        ggsave(png_str,
               get(str_glue("plt{x}_cov")),
               width = if (x == "GE") 18 else 9,
               height = 6)
    })
    
    list(G = pltG_cov,
         E = pltE_cov,
         P = pltP_cov,
         H = pltH_cov,
         GE = pltGE_cov)
}

# out <- plot_covariances("fb-final", 1, 1, "hpdi")
# out <- plot_covariances("fb-final-old", 1, 1, "hpdi")
