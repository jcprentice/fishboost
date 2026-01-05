{
    library(data.table)
    library(stringr)
    library(purrr)
    library(ggplot2)
    library(ggeasy)
    library(ggh4x)
    
    source("get_Tinfs.R")
}

compare_Tinfs <- function(dataset = "sim-test1", scen = 1, itn = 1) {

    data_dir <- str_glue("datasets/{dataset}/data")
    results_dir <- str_glue("datasets/{dataset}/results")
    gfx_dir <- str_glue("datasets/{dataset}/gfx")
    
    f <- str_glue("{results_dir}/scen-{scen}-{itn}.rds")
    if (!file.exists(f)) {
        stop("Can't find file ", f)
    }
    res <- readRDS(f)
    params <- res$params
    x <- res$popn[, .(id, sdp, trial, donor, group, Tinf)]
    
    y <- get_Tinfs(dataset, scen, itn) |>
        apply(1, quantile, na.rm = TRUE) |>
        t() |>
        as.data.table()
    y[, id := .I]
    
    xy <- merge(y, x[, .(id, sdp, trial, donor, Tinf)], by = "id")
    xy <- xy[sdp == "progeny" & donor == 0]
    setorder(xy, trial, Tinf, `50%`)
    xy[, id := .I]
    xy[, c("sdp", "donor") := NULL]
    xy2 <- melt(xy, id.vars = c("id", "trial"))
    xy2[, src := fifelse(variable == "Tinf", 2, 1)]
    
    subtitle_str <- with(params, str_c(dataset, "/", label, ": ",
                                       str_remove(description, ", convergence")))
    
    plt <- ggplot(xy) +
        geom_ribbon(aes(id, ymin = `0%`, ymax = `100%`),
                    fill = "lightgrey") +
        geom_ribbon(aes(id, ymin = `25%`, ymax = `75%`),
                    fill = "darkgrey") +
        geom_line(aes(id, `50%`), colour = "blue") +
        geom_line(aes(id, Tinf), colour = "red") +
        labs(x = "Individual",
             y = "Tinf",
             title = "True Tinf (red) vs inferred Tinf (blue)",
             subtitle = subtitle_str) +
        # scale_colour_manual(breaks = c("0%", "25%", "50%", "75%", "100%", "Tinf"),
        #                       values = c("blue", "lightblue", "white", "pink", "red", "black")) +
        facet_wrap(. ~ trial,
                   ncol = 2,
                   scales = "free",
                   labeller = labeller(
                       trial = c("1" = "Trial 1 Contacts",
                                 "2" = "Trial 2 Contacts"))) +
        theme_bw()
    plt
    ggsave(str_glue("datasets/{gfx_dir}/comparing_Tinf.png"),
           plt, width = 9, height = 6, unit = "in")
}

{
    dataset <- "sim-test3"
    scen <- 1
    itn <- 1
}

plt <- compare_Tinfs(dataset = dataset,
                     scen = scen, itn = itn)
