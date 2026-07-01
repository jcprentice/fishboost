{
    library(data.table)
    library(purrr)
    library(stringr)
    library(ggplot2)
    library(cowplot)

    source("utils.R")
    source("get_plot_matrix.R")
    source("rename_pars.R")
}

pars_boxplots <- function(dataset = "fb-test",
                          scens = 0,
                          st_str = "",
                          plt_shape = "traits",
                          output = "pdf") {
    if (FALSE) {
        dataset <- "fb-test"; scens <- 0
        st_str <- ""
        plt_shape <- "traits"; output <- "pdf"
    }

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
    } else {
        scens <- files |>
            str_extract("scen-(\\d+)-", 1) |>
            as.integer()
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
             variable.factor = FALSE)

    x[, scen := as.factor(str_c("s", scens[scen]))]

    pars <- unique(x$parameter)
    html_pars <- setNames(html_names(pars), pars)

    f <- list.files(str_glue("{res_dir}"), "scen", full.names = TRUE)[[1]]
    priors <- readRDS(f)$params$priors

    lims <- x[, .(min = min(value, na.rm = TRUE),
                  max = max(value, na.rm = TRUE)),
              parameter]

    plts <- map(pars, \(par) {
        if (FALSE) {
            i <- 1; par <- pars[[1]]
        }
        xp <- x[parameter == par]
        rng <- lims[parameter == par, .(min, max)] |>
            unlist() |> unname()

        ggplot(xp) +
            aes(x = scen, y = value) +
            geom_boxplot(fill = "tomato",
                         outliers = FALSE,
                         staplewidth = 1) +
            scale_y_continuous(limits = ~ range(.x, 0, rng)) +
            # coord_flip() +
            labs(x = "Scenario",
                 y = "Value",
                 title = html_pars[[par]]) +
            theme_classic() +
            theme(plot.title = element_markdown())
    }) |> setNames(pars)

    title_plt <- ggplot() +
        labs(title = str_glue("Dataset: '{dataset}'"),
             subtitle = st_str) +
        theme_classic() +
        theme(plot.title = element_text(size = 22),
              plot.subtitle = element_text(size = 16))

    plts$empty <- ggplot() + theme_classic()

    pmat <- get_plot_matrix(pars, plt_shape)

    plt <- plot_grid(title_plt,
                     plot_grid(plotlist = plts[pmat$plt_names],
                               nrow = pmat$nr,
                               ncol = pmat$nc,
                               align = "hv"),
                     ncol = 1,
                     rel_heights = c(0.5, pmat$nr))

    walk(output, \(op) {
        ggsave(str_glue("{gfx_dir}/{dataset}-all_boxplots.{op}"), plt,
               width = 3 * pmat$nc,
               height = 2 * (pmat$nr + 0.5))
    })

    plt
}

if (FALSE) {
    pars_boxplots("fb-test", c(1,2,7,8), "Testing BICI on FB data")
    # pars_boxplots("fb-final", 1:8, "Testing Unlinked vs Linked Traits and FEs")
    # pars_boxplots("fb-final2", 1:4, "Testing Unlinked vs Linked Traits and FEs")
    # pars_boxplots("fb-lp", 1:12, "Testing varying the LP")
    # pars_boxplots("fb-donors", 1:3, "Testing reclassifying Seeder fish")
    # pars_boxplots("fb-simple", 1:6, "Testing if we can fit the LP")
    # pars_boxplots("fb-simple-b", 1:6, "Testing if we can fit the LP with BICI")
}
