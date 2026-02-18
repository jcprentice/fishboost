{
    library(data.table)
    library(stringr)
    library(purrr)

    source("etc_to_km.R")
    source("generate_km_data_sire.R")
    source("generate_km_data_bici.R")
    source("run_bici_sim.R")
    source("plotting/plot_km_sire_trial.R")
    source("plotting/plot_km_donor_trial.R")
}

km_plots <- function(dataset = "fb-test",
                     scens = 0,
                     simulate_new_data = "bici",
                     opts = list(n_plots = 5,
                                 use_means = FALSE),
                     plotopts = c("drop_small_groups",
                                  "drop_donors",
                                  "use_sire_Tinfs")[0]) {

    if (FALSE) {
        dataset <- "sim-base-inf"; scens <- 1:18;  simulate_new_data <- "no"
        dataset <- "fb-qtest"; scens <- 0; simulate_new_data <- "bici"
        opts <- list(n_plots = 50, use_means = FALSE)
        plotopts <- c("drop_small_groups", "drop_donors", "extreme_sires")[1]
    }

    if (all(scens == 0)) {
        scens <- list.files(str_glue("datasets/{dataset}/results")) |>
            str_split_i("-", 2) |> unique() |> as.integer() |> sort()
    }

    # use_means <- FALSE

    # Use posterior means (pm) or sample from posterior (ps, default)
    um <- if (opts$use_means) "_pm" else "_ps"
    es <- if ("extreme_sires" %in% plotopts) "_es" else ""
    dd <- if ("drop_donors" %in% plotopts) "_dd" else ""

    f <- str_glue("datasets/{dataset}/meta/km_data{um}.rds")

    if (simulate_new_data != "no" && file.exists(f)) {
        message("- Removing previous km_data file")
        file.remove(f)
    }

    if (!file.exists(f)) {
        if (simulate_new_data == "no") {
            out <- etc_to_km(dataset, scens, opts, "post-sim")
            if (is.null(out)) {
                simulate_new_data <- "bici"
            }
        } else {
            simulate_new_data <- "bici"
        }

        if (simulate_new_data != "no") {
            message(str_glue("- No data found at '{f}', simulating new data ..."))
        }
    }

    if (simulate_new_data != "no") {
        message(str_glue("- Generating new data with '{simulate_new_data}'"))
        km_data <- if (simulate_new_data == "sire") {
            map(scens, ~ generate_km_data_sire(dataset, .x, opts))
        } else {
            map(scens, ~ generate_km_data_bici(dataset, .x, opts))
        }

        # Save the data we just created
        md <- str_glue("datasets/{dataset}/meta")
        if (!dir.exists(md)) {
            message("- mkdir ", md)
            dir.create(md)
        }
        saveRDS(km_data, file = f)
        message(str_glue("Saved KM data to '{f}'"))
    } else {
        km_data <- readRDS(f)
    }

    if (length(scens) != length(km_data)) {
        scens <- km_data |> map("params") |> map_int("scenario")
    }

    # Generate the plots
    plts_st <- map(km_data, plot_km_sire_trial,  plotopts)
    plts_dt <- map(km_data, plot_km_donor_trial, plotopts)

    # Grab the descriptions of the first available result for each scen
    descriptions <- map_chr(scens, \(i) {
        f <- list.files(str_glue("datasets/{dataset}/results"),
                        str_glue("scen-{i}-"),
                        full.names = TRUE)
        readRDS(f[[1]])$params$description
    })

    # Zip the plots with their descriptions
    plts_st <- map2(plts_st, descriptions, \(x, y) append(x, list(description = y)))
    plts_dt <- map2(plts_dt, descriptions, \(x, y) append(x, list(description = y)))

    km_dir_st <- str_glue("datasets/{dataset}/gfx/km_st")
    if (!dir.exists(km_dir_st)) {
        message("- mkdir ", km_dir_st)
        dir.create(km_dir_st, recursive = TRUE)
    }

    km_dir_dt <- str_glue("datasets/{dataset}/gfx/km_dt")
    if (!dir.exists(km_dir_dt)) {
        message("- mkdir ", km_dir_dt)
        dir.create(km_dir_dt, recursive = TRUE)
    }


    # Save the plots as PDF and PNG
    walk(seq_along(scens), \(i) {
        scen <- scens[[i]]
        plt_st <- plts_st[[i]]$plt
        st_str <- str_glue("{km_dir_st}/{dataset}-s{scen}-sire-trial{um}{es}{dd}")
        # ggsave(str_glue("{st_str}.pdf"), plt_st, width = 9, height = 6, unit = "in")
        ggsave(str_glue("{st_str}.png"), plt_st, width = 9, height = 6, unit = "in")
        message(str_glue("Saved plots to {st_str}"))

        plt_dt <- plts_dt[[i]]$plt
        dt_str <- str_glue("{km_dir_dt}/{dataset}-s{scen}-donor-trial{um}{es}{dd}")
        # ggsave(str_glue("{dt_str}.pdf"), plt_dt, width = 9, height = 6, unit = "in")
        ggsave(str_glue("{dt_str}.png"), plt_dt, width = 9, height = 6, unit = "in")
        message(str_glue("Saved plots to {dt_str}"))
    })
}

if (FALSE) {
    dataset <- "sim-base-inf"
    scens <- 0
    simulate_new_data <- "no"
    opts <- list(n_plots = 50, use_means = FALSE)
    plotopts <- c("drop_small_groups", "extreme_sires")[2]

    km_plots(dataset = dataset,
             scens = scens,
             simulate_new_data = simulate_new_data,
             opts = opts,
             plotopts = plotopts)
}

if (FALSE) {
    dataset <- "fb-qtest"
    scens <- 0
    simulate_new_data <- "bici"
    opts <- list(n_plots = 50, use_means = FALSE)
    plotopts <- c("drop_small_groups", "drop_donors", "extreme_sires")[1]

    km_plots(dataset = dataset,
             scens = scens,
             simulate_new_data = simulate_new_data,
             opts = opts,
             plotopts = plotopts)
}

if (FALSE) {
    dataset <- "fb-final"
    scens <- 1:8
    simulate_new_data <- "no"
    opts <- list(n_plots = 50, use_means = FALSE)
    plotopts <- "drop_small_groups"

    km_plots(dataset = dataset,
             scens = scens,
             simulate_new_data = simulate_new_data,
             opts = opts,
             plotopts = plotopts)
}

if (FALSE) {
    dataset <- "fb-final2"
    scens <- 1:4
    simulate_new_data <- "R"
    opts <- list(n_plots = 50, use_means = FALSE)
    plotopts <- "drop_small_groups"

    km_plots(dataset = dataset,
             scens = scens,
             simulate_new_data = simulate_new_data,
             opts = opts,
             plotopts = plotopts)
}

if (FALSE) {
    dataset <- "fb-lp"
    scens <- 1:12
    simulate_new_data <- "R"
    opts <- list(n_plots = 50, use_means = FALSE)
    plotopts <- "drop_small_groups"

    km_plots(dataset = dataset,
             scens = scens,
             simulate_new_data = simulate_new_data,
             opts = opts,
             plotopts = plotopts)
}

if (FALSE) {
    dataset <- "fb-lp2"
    scens <- 1:12
    simulate_new_data <- "R"
    opts <- list(n_plots = 50, use_means = FALSE)
    plotopts <- "drop_small_groups"

    km_plots(dataset = dataset,
             scens = scens,
             simulate_new_data = simulate_new_data,
             opts = opts,
             plotopts = plotopts)
}

if (FALSE) {
    dataset <- "fb-simple"
    scens <- 1:10
    simulate_new_data <- "R"
    opts <- list(n_plots = 50, use_means = FALSE)
    plotopts <- "drop_small_groups"

    km_plots(dataset = dataset,
             scens = scens,
             simulate_new_data = simulate_new_data,
             opts = opts,
             plotopts = plotopts)
}

if (FALSE) {
    dataset <- "fb-donors"
    scens <- 1:3
    simulate_new_data <- "R"
    opts <- list(n_plots = 50, use_means = FALSE)
    plotopts <- "drop_small_groups"

    km_plots(dataset = dataset,
             scens = scens,
             simulate_new_data = simulate_new_data,
             opts = opts,
             plotopts = plotopts)
}

if (FALSE) {
    dataset <- "sim-test1"
    scens <- 1:2
    simulate_new_data <- "R"
    opts <- list(n_plots = 50, use_means = FALSE)
    plotopts <- "drop_small_groups"

    km_plots(dataset = dataset,
             scens = scens,
             simulate_new_data = simulate_new_data,
             opts = opts,
             plotopts = plotopts)
}

if (FALSE) {
    opts <- list(n_plots = 50, use_means = FALSE)
    plotopts = c("drop_small_groups", "use_sire_Tinfs")

    km_plots(dataset = "fb-final",
             scens = 1,
             simulate_new_data = TRUE,
             opts = opts,
             plotopts = plotopts)

    km_plots(dataset = "sim-test2",
             scens = 1:4,
             simulate_new_data = FALSE,
             opts = opts,
             plotopts = plotsopts)
}


