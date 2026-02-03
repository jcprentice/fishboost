{
    library(data.table)
    library(purrr)
    library(stringr)
    library(lubridate)
    
    source("utils.R")
    source("rename_pars.R")
}

# Take the BICI state files and create an `etc.rds` file that contains all the
# parameters and traits in a data.table for easy access.

# bici_cmd: where do the output of BICI's simulations come from?
# "inf": look for all 'state_<n>.txt' files, and save to 'etc.rds'
# "sim": look for all 'state_<n>.txt' files, and save to 'etc-sim.rds'
# "post-sim": look for a single 'state.txt' file and save to 'etc-ps.rds'

flatten_bici_states <- function(dataset = "fb-test",
                                name = "scen-1-1",
                                bici_cmd = "inf") {
    # dataset <- "testing"; name <- "scen-1-1"; bici_cmd <- "inf"
    
    message(str_glue("Flattening BICI state files for '{dataset}/{name}' via '{bici_cmd}'"))
    
    base_dir <- str_glue("datasets/{dataset}")
    data_dir <- str_glue("{base_dir}/data")
    res_dir <- str_glue("{base_dir}/results")
    out_dir <- str_glue("{data_dir}/{name}-out")
    src_dir <- str_glue("{out_dir}/output-{bici_cmd}")
    
    rf <- str_glue("{res_dir}/{name}.rds")
    res <- readRDS(rf)
    base_popn <- res$popn[, .(id, sire, dam, sdp, trial, group, weight, donor)]
    params <- res$params
    
    files <- list.files(src_dir, "state", full.names = TRUE) |>
        str_sort(numeric = TRUE)
    
    if (is_empty(files)) {
        message("- No state files found!")
        return(invisible(NULL))
    }
    
    # Read first file to get positions of each state contained within
    lines <- readLines(files[[1]])
    
    # The top of the file contains the map of name:id. Transitions and IEs both
    # give `name`, where we want `id`, so just sub in `ids` where appropriate.
    idsect <- str_which(lines, "[{}]")
    ids <- lines[seq(idsect[[1]] + 1, idsect[[2]] - 1)] |>
        str_subset("#|Time|timepoint|index", negate = TRUE) |>
        str_squish() |>
        str_split_i(":", 2) |>
        as.integer()
    
    state_lines <- c(str_which(lines, "STATE"),
                     length(lines) + 1L)
    nstates <- length(state_lines) - 1L
    
    nfiles <- length(files)
    total <- nfiles * nstates
    
    message(str_glue("Found {nfiles} file{s} with {nstates} states, ",
                     "total = {total} states",
                     s = if (nfiles > 1) "s, each" else ","))
    
    # This will read all the lines between headers n and n+1
    get_section <- function(lines, headers, name) {
        n <- str_which(lines[headers], name)
        lines[seq(headers[[n]] + 1, headers[[n + 1]] - 1)]
    }
    
    parts <- c("parameters", "popn")
    
    t_start <- as.integer(now())
    
    out2 <- map(files, \(f) {
        # f <- files[[1]]
        
        lines <- readLines(f)
        
        out1 <- map(seq_len(nstates), \(i) {
            # i <- 1
            lines1 <- lines[seq(state_lines[[i]] + 2,
                                state_lines[[i + 1]] - 2)]
            
            headers <- c(str_which(lines1, "<.*>"),
                         length(lines1) + 1L)
            
            # Parameters
            p_lines <- get_section(lines1, headers, "PARAMETER") |>
                str_squish() |>
                str_subset(".+") |>
                str_subset("L\\^|Prior", negate = TRUE)
            
            p_headers <- c(str_which(p_lines, "\""),
                           length(p_lines) + 1L)
            
            parameters <- map(seq_len(length(p_headers) - 1), \(j) {
                # j <- 1
                phj  <- p_headers[[j]]
                phj1 <- p_headers[[j + 1]]
                
                par <- p_lines[phj] |> str_split_i("\"", 2)
                
                if (phj1 - phj > 1) {
                    p_tab <- fread(text = p_lines[seq(phj, phj1 - 1)] |>
                                       str_replace_all("\\.,", ","))
                    if (str_detect(par, "Omega")) {
                        GE <- str_split_i(par, "\\^", 2) |> str_1st() |> str_to_upper()
                        m <- t(p_tab)
                        traits <- rownames(m) |> str_1st()
                        tm <- expand.grid(traits, traits, "_", GE, "r_") |>
                            rev() |> apply(1, str_flatten) |>
                            matrix(length(traits))
                        diag(tm) <- str_replace(diag(tm), "r", "cov")
                        
                        p_tab <- data.table(parameter = c(diag(tm), tm[lower.tri(tm)]),
                                            Value = c(diag(m), m[lower.tri(m)]))
                    } else if ("g" %in% names(p_tab)) {
                        p_tab[, parameter := str_c(par, "_", g)]
                    } else if ("c" %in% names(p_tab)) {
                        p_tab[, parameter := str_c(par, "_", b, "," ,c)]
                    } else {
                        p_tab[, parameter := str_c(par, "_", b)]
                    }
                    p_tab[, .(parameter, value = Value)]
                } else {
                    value <- p_lines[[phj]] |>
                        str_split_i("\"", 3) |>
                        str_remove(",") |>
                        as.numeric()
                    data.table(parameter = par, value)
                }
            }) |>
                rbindlist()
            parameters[, parameter := rename_bici_pars(parameter)]
        
            
            # Transitions
            t_lines <- get_section(lines1, headers, "INDIVIDUALS") |>
                str_squish() |>
                str_subset(".+")
            
            get_event <- function(events, tr = "S->E") {
                if (str_starts(events[[1]], "unobserved")) {
                    return(NA)
                } else if (str_starts(tr, "S") && !str_starts(events[[1]], "S")) {
                    return(0)
                }
                ev <- str_subset(events, tr)
                if (length(ev) == 0) return(NA_real_)
                str_split_i(ev, ":", 2) |> as.numeric()
            }
            
            transitions <- fread(text = t_lines)[, .(id = ids, events)]
            iwalk(transitions$events, \(event, j) {
                # j <- 1; event <- transitions$events[[j]]
                events <- str_split_1(event, " ")
                switch(params$model_type,
                       "SEIDR" = {
                           set(transitions, j, c("Tinf", "Tinc", "Tsym", "Tdeath"),
                               list(Tinf = get_event(events, "S->E"),
                                    Tinc = get_event(events, "E->I"),
                                    Tsym = get_event(events, "I->D"),
                                    Tdeath = get_event(events, "D->R")))
                       }, "SIDR" = {
                           set(transitions, j, c("Tinf", "Tsym", "Tdeath"),
                               list(Tinf = get_event(events, "S->I"),
                                    Tsym = get_event(events, "I->D"),
                                    Tdeath = get_event(events, "D->R")))
                       }, "SEIR" = {
                           set(transitions, j, c("Tinf", "Tsym", "Tdeath"),
                               list(Tinf = get_event(events, "S->E"),
                                    Tinc = get_event(events, "E->I"),
                                    Tdeath = get_event(events, "I->R")))
                       }
                )
            })
            transitions[, events := NULL]
            popn <- merge(base_popn, transitions, by = "id", all = TRUE)

    
            if (params$use_traits %notin% c("none", "")) {
                i_lines <- get_section(lines1, headers, "INDIVIDUALS") |>
                    str_squish() |>
                    str_subset(".+")
                ies <- fread(text = i_lines)
                ies[, `:=`(index = ids, source = NULL, events = NULL)] |>
                    setnames(c("index",
                               "sg", "ig", "lg", "dg", "tg",
                               "se", "ie", "le", "de", "te"),
                             c("id",
                               "sus_g", "inf_g", "lat_g", "det_g", "tol_g",
                               "sus_e", "inf_e", "lat_e", "det_e", "tol_e"),
                             skip_absent = TRUE) |>
                    setorder(id)
                
                # inf returns LN dist'd IEs, need to rescale them
                if (ies[, min(.SD) > 0, .SDcols = -1]) {
                    ies[, names(.SD) := map(.SD, ~ log(.x) + var(log(.x)) / 2), .SDcols = -1]
                }
            } else {
                cols <- names(res$popn) |>
                    str_subset("id|_BV|_EV")
                cols2 <- cols |>
                    str_replace_all(c("susceptibility" = "sus",
                                      "infectivity" = "inf",
                                      "latency" = "lat",
                                      "detectability" = "det",
                                      "tolerance" = "tol",
                                      "BV" = "g", "EV" = "e"))
                ies <- res$popn[, ..cols] |>
                    setnames(cols, cols2)
            }
            
            popn <- merge(popn, ies, by = "id", all = TRUE)
            
            cat(".", file = stderr())
            
            mget(parts)
        })
        
        cat("\n", file = stderr())
        
        map(parts, ~ map(out1, .x) |> rbindlist(idcol = "state1")) |>
            setNames(parts)
    })
    t_end <- as.integer(now())
    message(str_glue("- Completed in {tt}",
                     tt = seconds_to_period(t_end - t_start)))
    
    etc <- map(parts, ~ map(out2, .x) |> rbindlist(idcol = "state2")) |>
        setNames(parts)
    
    map(etc, ~ {
        .x[, state := .GRP, .(state1, state2)]
        .x[, c("state1", "state2") := NULL]
        setcolorder(.x, "state")
    })
    
    
    f <- str_glue("{out_dir}/etc{x}.rds",
                  x = switch(bici_cmd,
                             "inf" = "_inf",
                             "sim" = "_sim",
                             "post-sim" = "_ps",
                             ""))
    message(str_glue("- Saving to '{f}'"))
    saveRDS(etc, f)
    
    invisible(etc)
}

# out <- flatten_bici_states("fb-test", "scen-1-1", "inf")
# out <- flatten_bici_states("testing", "scen-1-1", "sim")
# out <- flatten_bici_states("testing", "scen-1-1", "post-sim")
