{
    library(data.table)
    library(stringr)
    library(purrr)
    
    source("rename_pars.R")
}

dataset <- "fb-fes3"; scens <- 1:8

x <- map(scens, \(i) {
    f <- str_glue("datasets/{dataset}/results/scen-{i}-1.rds")
    PEs <- readRDS(f)$parameter_estimates[!str_starts(parameter, "Group effect")]
    PEs[ESS == "---", ESS := NA]
    PEs[, ESS := as.numeric(ESS)]
    PEs
}) |>
    rbindlist(idcol = "scen")

inf1 <- x[str_detect(parameter, "beta|cov_G|_i$") & scen %in% 1:4]
inf2 <- x[str_detect(parameter, "beta|cov_G|_i$") & scen %in% 5:8]
inf1[parameter == "beta"]
x[ESS < 200]


mtraces <- map(scens, \(i) {
    traces <- map(1:16, \(j) {
        # filter down to 100 samples, too many otherwise
        f <- str_glue("datasets/{dataset}/data/scen-{i}-1-out/trace_{j-1}.tsv")
        tmp <- fread(f)[round((1:100) * .N / 100)]
        tmp[, `:=`(h2_ss = cov_G_ss / (cov_G_ss + cov_E_ss),
                   h2_ii = cov_G_ii / (cov_G_ii + cov_E_ii),
                   h2_tt = cov_G_tt / (cov_G_tt + cov_E_tt))]
    }) |>
        rbindlist(idcol = "chain")
    
    traces[, str_subset(names(traces), "Group|L_|latent_period|Prior|Posterior|Number|log") := NULL]
    traces[, state := seq.int(0L, .N - 1L)]
    
    melt(traces, id.vars = c("chain", "state"))
}) |>
    rbindlist(idcol = "scen", fill = TRUE)

protocol <- readRDS(str_glue("param_sets/{dataset}.rds"))$protocol
descriptions <- protocol[, str_c("s", scenario, ": ", str_remove(description, ", convergence"))]

mtraces[, scenario_name := ordered(descriptions[scen])]

vars <- mtraces[, as.character(unique(variable))]
# FIXME this is ugly
vars <- intersect(param_order, vars)

walk(vars, \(v) {
    plt <- ggplot(mtraces[variable == v]) +
        geom_line(aes(x = state, y = value, colour = as.factor(chain))) +
        labs(title = str_glue("Parameter: {v}"), colour = "Chain") +
        facet_wrap(. ~ scenario_name, nrow = 2) +
        theme(legend.position = "none")
    
    ggsave(str_glue("datasets/{dataset}/gfx/{dataset}-trace-{v}.png"),
           plt, width = 12, height = 9)
})
