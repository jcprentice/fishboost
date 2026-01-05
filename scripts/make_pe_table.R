{
    library(data.table)
    library(purrr)
    library(stringr)
    library(HDInterval)
    library(xtable)
    
    source("rename_pars.R")
}

# dataset <- "fb-fes3"; scens <- 1:12
# dataset <- "fb-final-old"; scens <- 1:8
dataset <- "fb-final"; scens <- 1:4
# dataset <- "fb-simple"; scens <- 1:10
# dataset <- "sim-test"; scens <- 1:4
# dataset <- "sim-final"; scens <- 1:2

scens <- setNames(scens, str_c("s", scens))

x <- map(scens, possibly(\(i) {
    # i <- 1
    f <- str_glue("datasets/{dataset}/data/scen-{i}-1-out/trace_combine.tsv")
    if (!file.exists(f)) return(NULL)
        
    tr <- fread(f)
    tr[, str_subset(names(tr), "state|Group effect") := NULL]
    if (any(str_detect(names(tr), "cov_G"))) {
        tr[, `:=`(h2_ss = cov_G_ss / (cov_G_ss + cov_E_ss),
                  h2_ii = cov_G_ii / (cov_G_ii + cov_E_ii),
                  h2_tt = cov_G_tt / (cov_G_tt + cov_E_tt))]
    }
    
    mu <- tr[, map(.SD, mean)]
    hpdi <- hdi(tr, credMass = 0.95)
    
    pe <- data.table(parameter = names(mu),
                     mean = unlist(mu),
                     min95 = hpdi["lower", ],
                     max95 = hpdi["upper", ])
    
    rows <- pe[, str_subset(parameter, "h2", negate = TRUE)]
    
    f <- str_glue("datasets/{dataset}/results/scen-{i}-1.rds")
    if (file.exists(f)) {
        PEs <- readRDS(f)$parameter_estimates
        pe[parameter %in% rows, `:=`(ESS = PEs[parameter %in% rows, ESS],
                                     GR = PEs[parameter %in% rows, GR])]
    } else {
        pe[, `:=`(ESS = NA_real_, GR = NA_real_)]
    }
    
    if (is.character(pe$ESS)) stop("Need to rebuild posteriors.csv and fix results")
    
    setnames(pe, c("parameter", "mean", "min95", "max95", "ESS", "GR"))
    # pe[parameter == "latent_period", `:=`(ESS = NA, GR = NA)]
    # pe[, ESS := as.numeric(ESS)]
    pe
})) |>
    rbindlist(idcol = "i") |>
    dcast(parameter ~ i, value.var = c("mean", "min95", "max95", "ESS", "GR"))

# Fix names if they don't line up with scens
tmp <- expand.grid(scen = scens, var = c("mean", "min95", "max95", "ESS", "GR"))[, c("var", "scen")] |>
    apply(1, str_flatten, "_") |>
    str_remove_all(" ")
setnames(x, c("parameter", tmp))

# The order of the colnames is [parameter, mean1, mean2, ... , GR7, GR8]. The
# simplest way I can see to fix that is to put the colnames into a matrix and
# transpose it.

ncols <- sum(str_detect(names(x), "mean"))
setcolorder(x, c("parameter", t(matrix(names(x)[-1], nrow = ncols))))


# Rename_pars.R gives the values needed here
x[, row_order := match(parameter, param_order)]
setorder(x, row_order)
x[, row_order := NULL]
x[, parameter := rename_pars(parameter)]


cols <- names(x)[-1]
x[, (cols) := round(.SD, 3), .SDcols = cols]

# Collapse min95 and max95 into single string
walk(scens, possibly(\(i) {
    min_i <- str_c("min95_", i)
    max_i <- str_c("max95_", i)
    hpdi95_i <- str_c("hpdi95_", i)
    x[, (hpdi95_i) := str_c("(", get(min_i), ", ", get(max_i), ")")]
    setcolorder(x, hpdi95_i, after = str_c("mean_", i))
    x[, c(min_i, max_i) := NULL]
}))

setnames(x, str_replace_all(names(x), "_", "_s"))

pe_name <- str_glue("dataset/meta/{dataset}-parameter_estimates.txt")
fwrite(x, pe_name, sep = "|")

# setnames(x, c("parameter", rep(cols, 3)))

if (FALSE) {
    hlines <- c(-1, 0,
                head(which(str_detect(x$parameter, "Var")), 1) - 1,
                head(which(str_detect(x$parameter, "Trial|Donor|TxD|Weight")), 1) - 1,
                nrow(x))
    
    print(xtable(x, align = c("l", "l", rep("r", ncol(x) - 1))),
          include.rownames = FALSE,
          booktabs = TRUE,
          size = "small",
          hline.after = hlines)
}

