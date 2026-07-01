{
    library(data.table)
    library(stringr)
    library(purrr)
    library(ggplot2)
    library(ggtext)
    library(ggeasy)
    library(RColorBrewer)
    library(cowplot)

    source("utils.R")
    source("rename_pars.R")
}

plot_chains <- function(dataset = "fb-test",
                        scen = 1, rep = 1,
                        plt_shape = "traits") {
    if (FALSE) {
        dataset <- "sim-test-inf1"; scen <- 1; rep <- 1
        dataset <- "fb-test"; scen <- 7; rep <- 1
        plt_shape <- "traits"
    }

    message(str_glue("plotting chain for '{dataset}' / s{scen}-{rep}"))

    {
        base_dir <- str_glue("datasets/{dataset}")
        data_dir <- str_glue("{base_dir}/data")
        res_dir <- str_glue("{base_dir}/results")
        chains_dir <- str_glue("{base_dir}/gfx/chains")
    }

    if (!dir.exists(chains_dir)) {
        message("- mkdir ", chains_dir)
        dir.create(chains_dir, recursive = TRUE)
    }

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
        priors <- priors[parameter %in% pars,
                         .(parameter, min = val1, max = val2, true_val)]
        description <- res_file$params$description
    } else {
        priors <- ltraces[, .(min = min(y), max = max(y),
                              true_val = NA_real_),
                          parameter]
        description <- ""
    }

    nchains <- ltraces$chain |> levels() |> length()
    thin <- 1L
    plts <- map(pars, \(par) {
        if (FALSE) {
            i <- 1
            par <- pars[[i]]
        }
        X <- ltraces[parameter == par][seq(1, .N, thin)]
        ytrue <- if (str_detect(dataset, "sim")) {
            priors[parameter == par, true_val]
        }

        yrng <- ltraces[parameter == par, range(ytrue, y)]
        yrng <- yrng + c(-0.1, 0.1) * abs(yrng)

        p_true <- if (str_detect(dataset, "sim")) {
            geom_hline(yintercept = ytrue,
                       linetype = "dashed",
                       colour = "green")
        }

        lp2 <- 1 / ggplot2::.pt

        ggplot(X, aes(x, y, group = par, colour = chain)) +
            geom_point(size = 0.8 * lp2) +
            p_true +
            scale_colour_manual(values = rep(brewer.pal(4, "Set1"),
                                             length = nchains)) +
            scale_y_continuous(limits = ~ range(.x, 0, yrng)) +
            labs(x = NULL,
                 y = "value",
                 title = html_names(par)) +
            theme_classic() +
            theme(plot.title = element_markdown(size = 12),
                  plot.subtitle = element_markdown(size = 10),
                  legend.position = "none",
                  axis.text.x = element_blank())
    }) |>
        setNames(pars)

    title_plt <- ggplot() +
        labs(title = str_glue("{dataset} / s{scen}-{rep}"),
             subtitle = description) +
        theme_classic()
        # theme_classic() +
        # theme(plot.title = element_text(size = 12),
        #       plot.subtitle = element_text(size = 10))

    plts$empty <- ggplot() + theme_classic()

    pmat <- get_plot_matrix(pars, plt_shape)

    plt <- plot_grid(title_plt,
                     plot_grid(plotlist = plts[pmat$plt_names],
                               nrow = pmat$nr,
                               ncol = pmat$nc,
                               align = "v"),
                     ncol = 1, rel_heights = c(0.5, pmat$nr))

    # PDFs are huge here
    plt_str <- str_glue("{chains_dir}/{dataset}-s{scen}-{rep}-chains.png")
    ggsave(plt_str, plt,
           width = 4 * pmat$nc,
           height = 2 * (pmat$nr + 0.5))
    plt
}

if (FALSE) {
    dataset <- "sim-test-inf1"
    files <- str_glue("datasets/{dataset}/results") |>
        list.files() |> str_sort(numeric = TRUE)
    scens <- files |> str_extract("scen-(\\d+)-", 1) |> as.integer()
    reps <- files |> str_extract("scen-\\d+-(\\d+).rds", 1) |> as.integer()

    walk2(scens, reps, \(scen, rep) plot_chains(dataset, scen, rep))
}
