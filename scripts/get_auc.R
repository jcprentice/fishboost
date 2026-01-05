{
    library(data.table)
    library(purrr)
    library(stringr)
    library(ggplot2)
    library(ggeasy)
    library(ggh4x)
    
    # source("scripts/km_plots.R")
}

get_auc <- function(dataset = "fb-final",
                    gps = c("trial", "donor"),
                    opts = list(use_means = FALSE,
                                apply_gammas = 0)) {
    
    um <- if (opts$use_means) "_pm" else "_ps" # sample from prior
    ag <- if (opts$apply_gammas > 0) str_glue("_gam{opts$apply_gammas}") else ""
    
    f <- str_glue("datasets/{dataset}/meta/km_data{um}{ag}.rds")
    if (!file.exists(f)) {
        message(str_glue("'{f}' not found, generating now"))
        km_plots(dataset = dataset,
                 scens = 0,
                 simulate_new_data = TRUE,
                 opts = list(n_plots = 50,
                             use_means = FALSE,
                             apply_gammas = 0),
                 plotopts = list("drop_small_groups"))
    }
    km_data <- readRDS(f)
    
    map(km_data, \(x) {
        # x <- km_data[[1]]
        
        # Make a copy of the data and add in rows representing t=0, sort times
        x1 <- x$data[, .(Tsym = c(0, sort(Tsym, na.last = TRUE)),
                         RP   = c(0, sort(RP,   na.last = TRUE)),
                         src = src[c(1, 1:.N)]),
                     mget(c("id", gps))] |>
            setorderv(c("id", gps))
        
        # Create a column to group by, ensure this is wide enough to handle all ids
        x1[, id1 := format(id)]
        x1[, gp := apply(.SD, 1, str_flatten, " "), .SDcols = c("id1", gps)]
        x1[, id1 := NULL]
        
        
        # Survival curves by id and sire
        x1[, survival := seq(1, 0, length.out = .N), gp]
        
        # Melt so we can use facet_wrap
        x2 <- x1 |> melt(measure.vars = c("Tsym", "RP"))
        
        # Censor the data as with the FB data
        tmax <- c(104, 160) # data[src == "fb", ceiling(max(Tdeath, na.rm = TRUE)), trial][, V1]
        # x2[, value := fifelse(value < tmax[trial], value, tmax[trial], tmax[trial])]
        
        auc <- x2[, {
            v1 <- seq(0, tmax[trial])
            s1 <- approxfun(value, survival, ties = list("ordered", max))(v1)
            s1[is.na(s1)] <- 0
            list(value = v1, survival = s1)
        }, c("id", gps, "variable")]
        
        # This relies on FB being last
        auc[, auc_diff := survival - survival[[.N]], value]
        auc2 <- auc[, .(auc_diff = sum(abs(auc_diff))), c("id", gps, "variable")]
        
        auc2
    }) |>
        rbindlist(idcol = "scen")
}


get_auc_donor <- function(dataset = "fb-final",
                          opts = list(use_means = FALSE,
                                      apply_gammas = 0)) {
    
    um <- if (opts$use_means) "_pm" else "_ps" # sample from prior
    ag <- if (opts$apply_gammas > 0) str_glue("_gam{opts$apply_gammas}") else ""
    
    f <- str_glue("datasets/{dataset}/meta/km_data{um}{ag}.rds")
    if (!file.exists(f)) {
        message(str_glue("'{f}' not found, generating now"))
        km_plots(dataset = dataset,
                 scens = 0,
                 simulate_new_data = TRUE,
                 opts = list(n_plots = 50,
                             use_means = FALSE,
                             apply_gammas = 0),
                 plotopts = list("drop_small_groups"))
    }
    km_data <- readRDS(f)
    
    imap(km_data, \(x, scen) {
        data <- x$data
        
        # Make a copy of the data and add in rows representing t=0, sort times
        x1 <- data[, .(Tsym = c(0, sort(Tsym, na.last = TRUE)),
                       RP   = c(0, sort(RP,   na.last = TRUE)),
                       src = src[c(1, 1:.N)]),
                   .(id, donor, trial)]
        
        setkey(x1, id, donor, trial)
        
        # Create a column to group by, ensure this is wide enough to handle all ids
        x1[, gp := str_c(format(id), " ", donor, " ", trial)]
        
        # Survival curves by id and sire
        x1[, survival := seq(1, 0, length.out = .N), gp]
        
        # Melt so we can use facet_wrap
        x2 <- x1 |> melt(measure.vars = c("Tsym", "RP"))
        
        # Censor the data as with the FB data
        tmax <- c(104, 160)
        # x2[, value := fifelse(value < tmax[trial], value, tmax[trial], tmax[trial])]
        
        auc <- x2[, {
            v1 <- seq(0, tmax[trial])
            s1 <- approxfun(value, survival, ties = list("ordered", max))(v1)
            s1[is.na(s1)] <- 0
            list(value = v1, survival = s1)
        }, .(id, trial, donor, variable)]
        
        # This relies on FB being last
        auc[, auc_diff := survival - survival[[.N]], value]
        auc2 <- auc[, .(auc_diff = sum(abs(auc_diff))), .(id, trial, donor, variable)]
        
        auc2
    }) |>
        rbindlist(idcol = "scen")
}


get_auc_sire <- function(dataset = "fb-final",
                         opts = list(use_means = FALSE,
                                     apply_gammas = 0)) {
    
    um <- if (opts$use_means) "_pm" else "_ps" # sample from prior
    ag <- if (opts$apply_gammas > 0) str_glue("_gam{opts$apply_gammas}") else ""
    
    f <- str_glue("datasets/{dataset}/meta/km_data{um}{ag}.rds")
    if (!file.exists(f)) {
        message(str_glue("'{f}' not found"))
    }
    km_data <- readRDS(f)
    
    imap(km_data, \(x, scen) {
        data <- x$data
        
        # Make a copy of the data and add in rows representing t=0, sort times
        x1 <- data[, .(Tsym = c(0, sort(Tsym, na.last = TRUE)),
                       RP   = c(0, sort(RP,   na.last = TRUE)),
                       src = src[c(1, 1:.N)]),
                   .(id, sire, trial)]
        
        setkey(x1, id, sire, trial)
        
        # Create a column to group by, ensure this is wide enough to handle all ids
        x1[, gp := str_c(format(id), " ", sire, " ", trial)]
        
        # Survival curves by id and sire
        x1[, survival := seq(1, 0, length.out = .N), gp]
        
        # Melt so we can use facet_wrap
        x2 <- melt(x1,
                   measure.vars = c("Tsym", "RP"))
        
        # Censor the data as with the FB data
        tmax <- c(104, 160)
        x2[, value := fifelse(value < tmax[trial],
                              value, tmax[trial], tmax[trial])]
        
        # auc <- x2[, .(auc = sum(diff(value) * survival[-.N])),
        #                .(trial, sire, src, gp, variable)]
        # auc[, auc_diff := auc - auc[[.N]], .(variable, trial, sire)]
        
        auc <- x2[, {
            v1 <- seq(0, tmax[trial[[1]]])
            s1 <- approxfun(value, survival, ties = list("ordered", max))(v1)
            s1[is.na(s1)] <- 0
            list(value = v1, survival = s1)
        }, .(id, sire, trial, src, variable)]
        
        # This relies on FB being last
        auc[, auc_diff := survival - survival[[.N]], value]
        auc2 <- auc[, .(auc_diff = sum(abs(auc_diff))), .(id, sire, trial, variable)]
        
        auc2
    }) |>
        rbindlist(idcol = "scen")
}

dataset <- "fb-final"
opts = list(use_means = FALSE,
            apply_gammas = 0)

auc_all   <- get_auc(dataset, c("trial", "donor"), opts)
auc_donor <- get_auc_donor(dataset, opts)
auc_sire  <- get_auc_sire(dataset, opts)

scen_levels <- str_c("s", sort(unique(auc_donor$scen)))

auc_donor[, `:=`(scen = factor(str_c("s", scen), scen_levels),
                 td = str_c(trial, "_", donor))]
auc_sire[, scen := factor(str_c("s", scen), scen_levels)]
auc_all[, `:=`(scen = factor(str_c("s", scen), scen_levels),
               td = str_c(trial, "_", donor))]


saveRDS(mget(c("auc_all", "auc_donor", "auc_sire")),
        str_glue("datasets/{dataset}/meta/auc_data.rds"))

auc_donor_order <- auc_donor[, .(diff = mean(auc_diff)), scen]
auc_donor_order <- auc_donor[scen %in% str_c("s", 1:6), .(diff = mean(auc_diff)), scen]
auc_donor_order[, order := order(diff)]
auc_donor_order[, diff := round(diff, 1)]
auc_donor_order

plt_donor <- auc_ |>
    ggplot(aes(scen, auc_diff, fill = variable)) +
    geom_boxplot(outliers = FALSE,
                 staplewidth = 0.5,
                 show.legend = FALSE) +
    # geom_jitter(size = 0.02) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(x = "Scenario",
         y = "Difference in AUC from FB",
         title = str_glue("Dataset: '{dataset}', by donor / trial")) +
    expand_limits(y = 0) +
    theme_bw() +
    facet_grid2(rows = vars(td),
                cols = vars(variable),
                scales = "free",
                independent = "y",
                labeller = labeller(
                    td = c("1_0" = "Trial 1 Contact",
                           "1_1" = "Trial 1 Seeder",
                           "2_0" = "Trial 2 Contact",
                           "2_1" = "Trial 2 Seeder"),
                    variable = c("Tsym" = "Time to symptoms",
                                 "RP" = "Time from symptoms to death")))
plt_donor

ggsave(str_glue("datasets/{dataset}/gfx/auc_donor.png"),
       plt_donor, width = 8, height = 5)

plt_sire <- auc_sire |>
    ggplot(aes(scen, auc_diff, fill = variable)) +
    geom_boxplot(outliers = FALSE,
                 staplewidth = 0.5,
                 show.legend = FALSE) +
    # geom_jitter(size = 0.02) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(x = "Scenario",
         y = "Abs difference in AUC",
         title = str_glue("Dataset: '{dataset}', by sire / trial")) +
    expand_limits(y = 0) +
    theme_bw() +
    facet_grid2(rows = vars(trial),
                cols = vars(variable),
                scales = "free",
                independent = "y",
                labeller = labeller(
                    trial = c("1" = "Trial 1",
                              "2" = "Trial 2"),
                    variable = c("Tsym" = "Time to symptoms",
                                 "RP" = "Time from symptoms to death")))
plt_sire

ggsave(str_glue("datasets/{dataset}/gfx/auc_sire.png"),
       plt_sire, width = 8, height = 5)

