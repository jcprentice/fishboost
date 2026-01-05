{
    library(data.table)
    library(purrr)
    library(coda)
    library(xtable)
    library(stringr)
    
    source("rename_pars.R")
}


gr_diagnostics <- function(dataset = "fb-final", scen = 1, itn = 0) {
    
    # dataset <- "fb-final"; scen <- 1; itn <- 0
    
    data_dir <- str_glue("datasets/{dataset}/data")
    res_dir  <- str_glue("datasets/{dataset}/results")
    
    trace_files <- list.files(path = data_dir,
                              pattern = "trace",
                              recursive = TRUE,
                              full.names = TRUE) |>
        # Sieve for the right scenario
        str_subset(str_glue("scen-{scen}-")) |>
        # Drop combine files for SIRE 2.2
        str_subset("combine", negate = TRUE)
    
    # If doing coverage, then specify which itn to keep, or keep them all
    if (itn == 0) {
        trace_files <- trace_files |>
            str_subset(str_glue("-{itn}-out")) |>
            str_sort(numeric = TRUE)
    }
    
    f <- str_glue("{res_dir}/scen-{scen}-{itn}.rds")
    # f <- list.files(res_dir, pattern = str_glue("scen-{scen}-"), full.names = TRUE)
    
    # if (!file.exists(res[1])) {
    #     return(NULL)
    # }
    res <- readRDS(f[[1]])
    PEs <- res$parameter_estimates
    params <- res$params
    
    fixed_pars <- params$priors[type == "Fixed", parameter]
    if (is_empty(fixed_pars)) fixed_pars <- "don't_match"
    
    # Select columns to keep
    cols <- PEs$parameter |>
        str_subset("Group effect", negate = TRUE) |>
        # Remove columns with zero variance (e.g. "latent_period")
        str_subset(fixed_pars, negate = TRUE) # |>
        # str_subset("L_|state|Prior|Posterior|Number infected|log(phi)", negate = TRUE)
    
    
    thin <- with(params, if (exists("thinto")) nsample / thinto else thin)
    burnin <- with(params, if (exists("burnprop")) ceiling(nsample * burnprop) else burnin)
    bstart <- ceiling(burnin / thin + 1L)
    
    x <- map(trace_files, fread)
    # In case we didn't finish, all chains have to be the same length
    max_rows <- map_int(x, nrow) |> min()
    l <- map(x, ~ mcmc(.x[seq_len(max_rows), ..cols], start = bstart, thin = thin)) |>
        as.mcmc.list()
    y <- map(x, ~ as.data.table(.x[bstart:.N, ..cols])) |>
        rbindlist()
    rm(x)
    
    list(ESS = effectiveSize(y),
         GRD = gelman.diag(l, autoburnin = FALSE))
}

# dataset <- "fb-final-old"; scens <- 1:8
dataset <- "fb-final"; scens <- 1:8
# dataset <- "fb-final2"; scens <- 1:4
# dataset <- "sim-test1"; scens <- 1:9
# dataset <- "fb-simple"; scens <- 1:10
# dataset <- "fb-lp"; scens <- 1:12
dataset <- "fb-donors"; scens <- 1:3

itn <- 1 # which iteration if doing coverage

# Add in scenario names
scens <-  set_names(scens, str_c("s", scens))

# Perform GR diagnostics on each scenario
out <- map(scens, possibly(~ gr_diagnostics(dataset, .x, itn)))

# Extract the MPSRF for each scenario
gd <- map_dbl(out, \(x) {
    if (is.null(x)) return(NA)
    x$GRD$mpsrf |> round(3)
})
gd



# Extract the GR95 and ESS values
xl <- map(out, \(x) {
    if (is.null(x)) return()
    y <- as.data.table(x$GRD$psrf, keep.rownames = "Parameter") |>
        setnames(c("Parameter", "pe", "GR95"))
    y[, pe := NULL]
    y[, ESS := format(ceiling(x$ESS), big.mark = ",")]
    y[, Parameter := rename_pars(Parameter)]
    rbind(y, data.table(Parameter = "MVPSF", GR95 = round(x$GRD$mpsrf, 2)), fill = TRUE)
}) |> discard(is.null)


# Create a LaTeX table with a subset of scenarios
{
    # I'd like to do a reduce / merge, but I can't quite figure out how to avoid
    # duplicate col names, but this seems to work even better.
    xl_tab <- xl |>
        rbindlist(idcol = "scen") |>
        dcast(... ~ scen, value.var = c("ESS", "GR95"), drop = FALSE)
    
    # Restore the preferred parameter order (see 'rename_pars.R')
    xl_tab[, ids := match(Parameter, full_param_order)]
    setorder(xl_tab, ids)
    xl_tab[, ids := NULL]
    
    # Reorder by scenario, rather than variable
    cols <- expand.grid(c("ESS", "GR95"), str_c("s", scens)) |>
        apply(1, str_flatten, "_") |>
        intersect(names(xl_tab))
    setcolorder(xl_tab, c("Parameter", cols))
    
    align_str <- str_flatten(c("ll", rep("rr", length(cols) / 2)))
    
    xl_ltab <- print(xtable(xl_tab, align = align_str),
                     include.rownames = FALSE,
                     booktabs = TRUE)
    
    print(xtable(xl_tab, align = align_str),
          include.rownames = FALSE,
          booktabs = TRUE,
          file = "xl_tab.tex")
}

xl_tab[, map(.SD, min, na.rm = TRUE), .SDcols = str_subset(names(xl_tab), "ESS")]
xl_tab[, map(.SD, max, na.rm = TRUE), .SDcols = str_subset(names(xl_tab), "GR")]
gd

results <- c("dataset", "scens", "out", "gd", "xl", "xl_tab", "xl_ltab") |> keep(exists)

saveRDS(mget(results),
        file = str_glue("{dataset}/meta/gr.rds"))

