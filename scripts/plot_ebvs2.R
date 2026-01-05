{
    library(data.table)
    library(purrr)
    library(stringr)
    library(ggplot2)
    library(cowplot)
}

plot_ebvs2 <- function(dataset = "fb-qtest", scen = 1, rep = 1) {
    if (FALSE) {
        dataset <- "fb-qtest"; scen <- 1; rep <- 1
    }
    
    message(str_glue("Plotting ebvs2 for {dataset}/s{scen}-{rep}"))
    
    {
        base_dir <- str_glue("datasets/{dataset}")
        data_dir <- str_glue("{base_dir}/data")
        res_dir  <- str_glue("{base_dir}/results")
        gfx_dir  <- str_glue("{base_dir}/gfx")
        ebvs_dir <- str_glue("{gfx_dir}/ebvs2")
        
        if (!dir.exists(ebvs_dir)) dir.create(ebvs_dir)
    }
    
    rf <- str_glue("{res_dir}/scen-{scen}-{rep}.rds")
    desc <- readRDS(rf)$params$description |>
        str_remove_all(", convergence|, GRM \\w*") |>
        str_replace_all(c("inf_model 1" = "inf: I = D",
                          "inf_model 2" = "inf: I = 0.1*D",
                          "inf_model 3" = "inf: Don = 0.1*Rec",
                          "inf_model 4" = "inf: Don = r*Rec"))
    
    x <- readRDS(str_glue("{data_dir}/scen-{scen}-{rep}-out/etc_inf.rds"))
    
    BVs <- x$popn[sdp == "sire" | (sdp == "progeny" & !is.na(sire)),
                  c(list(sdp = as.character(first(sdp)),
                         trial = first(trial),
                         donor = first(donor)),
                    map(.SD, median)),
                  id,
                  .SDcols = c("sus_g", "inf_g", "tol_g")]
    BVs[sdp == "sire", grp := "Sire"]
    # BVs[sdp == "dam", grp := "Dam"]
    BVs[sdp == "progeny", grp := str_c("Trial ", trial, " ", fifelse(donor == 1, "Don", "Rec"))]
    setkey(BVs, NULL)
    
    groups <- BVs$grp |> unique() |> str_sort()
    values <- c("green", "red", "pink", "blue", "lightblue")
    labels <- groups
    
    muS <- median(BVs$sus_g)
    muI <- median(BVs$inf_g)
    muT <- median(BVs$tol_g)
    
    pSIa <-  ggplot(BVs, aes(x = sus_g, y = inf_g, colour = grp)) +
        geom_point() +
        geom_vline(xintercept = muS, linetype = "dashed") +
        geom_hline(yintercept = muI, linetype = "dashed") +
        scale_colour_manual("Group", breaks = groups, values = values, labels = labels) +
        labs(x = "Susceptibility",
             y = "Infectivity") +
        theme_bw()
    pSI <- pSIa + theme(legend.position = "none")
    pSIa
    legend <- get_legend(pSIa)
    
    pST <- ggplot(BVs, aes(x = sus_g, y = tol_g, colour = grp)) +
        geom_point() +
        geom_vline(xintercept = muS, linetype = "dashed") +
        geom_hline(yintercept = muT, linetype = "dashed") +
        scale_colour_manual("Group", breaks = groups, values = values, labels = labels) +
        labs(colour = "Group",
             x = "Susceptibility",
             y = "Tolerance") +
        theme_bw() +
        theme(legend.position = "none")
    pST
    
    pIT <- ggplot(BVs, aes(x = inf_g, y = tol_g, colour = grp)) +
        geom_point() +
        geom_vline(xintercept = muI, linetype = "dashed") +
        geom_hline(yintercept = muT, linetype = "dashed") +
        scale_colour_manual("Group", breaks = groups, values = values, labels = labels) +
        labs(colour = "Group",
             x = "Infectivity",
             y = "Tolerance") +
        theme_bw() +
        theme(legend.position = "none")
    pIT
    
    title_plt <- ggplot() +
        labs(title = str_glue("{dataset} / s{scen}-{rep}"),
             subtitle = desc) +
        theme_classic() +
        theme(plot.title = element_text(size = 22),
              plot.subtitle = element_text(size = 16))
    
    
    pSIT <- plot_grid(title_plt,
                      plot_grid(pSI, pST, pIT, legend),
                      ncol = 1,
                      rel_heights = c(0.06, 1))
    pSIT
    
    plot_str <- str_glue("{ebvs_dir}/{dataset}-s{scen}-{rep}-ebvs2.pdf")
    ggsave(plot_str, pSIT, width = 12, height = 12)
    
    pSIT
}

if (FALSE) {
    dataset <- "fb-qtest"
    scens <- list.files(str_glue("datasets/{dataset}/results")) |>
        str_split_i("-", 2) |>
        as.integer() |>
        sort()
    
    walk(scens, ~ plot_ebvs2("fb-qtest", .x, 1))
}
