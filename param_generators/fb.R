{
    library(data.table)
    library(purrr)
    library(stringr)

    source("utils.R")
}

# Overview ----

# Are we testing convergence or coverage?
# (Coverage only makes sense with simulated data.)
n <- 1
goal <- c("convergence", "coverage")[[n]]
message("Goal = ", goal)

dataset <- "fb"

# Variable parameters ----
protocol <- rbind(
    data.table(d = "FB_12, SEIDR, Traits SIT, FEs ILDT, GE 0.3"),
    data.table(d = "FB_1,  SEIDR, Traits SIT, FEs ILDT, GE 0.3"),
    data.table(d = "FB_2,  SEIDR, Traits SIT, FEs ILDT, GE 0.3"),
    data.table(d = "FB_1, GE -1",
               trial_fe = ""),
    
    fill = TRUE
)

protocol[, setup := d |> str_split_i(", ", 1) |> str_to_lower()]

protocol[, group_effect := get_part(d, "GE") |> as.numeric(), .I]

# Common options ----
source("param_generators/common2.R")

common <- list(use_traits = "sit",
               donor_fe = "ildt") |>
    safe_merge(common2)

# Labels
protocol[, label := str_c("s", 1:.N)]

# Append "coverage" or "convergence" to description
protocol[, d := str_c(d, ", ", goal) |> str_squish()] |>
    setnames("d", "description")

## Add replicates ----
n_replicates <- if (goal == "convergence") 1 else 20
protocol[, scenario := .I]
protocol <- protocol[rep(1:.N, each = n_replicates)]
protocol[, replicate := 1:.N, scenario]
protocol[, dataset := dataset]
protocol[, name := str_c("scen-", scenario, "-", replicate)]

# Prefer to have these columns in this order at the start
setcolorder(protocol, intersect(cols, names(protocol)))

message(str_glue("Protocol file '{dataset}' has:",
                 "- {nrow(protocol)} scenarios",
                 "- each with {ncol(protocol)} parameters",
                 "- and a further {length(common)} common parameters",
                 .sep = "\n"))

# Save to file ----
saveRDS(list(protocol = protocol,
             common = common),
        file = str_glue("param_sets/{dataset}.rds"))

