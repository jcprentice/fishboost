source("models/model_SIS.R")
source("models/model_SIR.R")
source("models/model_SEIR.R")
source("models/model_SEIR_res.R")
source("models/model_SIDR.R")
source("models/model_SEIDR.R")

simulate_epidemic <- function(popn, params) {
    # this just turns "simulate_epidemic(...) into "model_SIR(...)"
    simulate_fn <- get(str_c("model_", params$model_type))
    simulate_fn(popn, params)
}

