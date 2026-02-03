# source("models/model_SI.R")
source("models/model_SIS.R")
source("models/model_SIR.R")
source("models/model_SEIR.R")
source("models/model_SEIR_res.R")
source("models/model_SIDR.R")
source("models/model_SEIDR.R")

simulate_epidemic <- function(popn, params) {
    switch(params$model_type,
           "SI"  = model_SI(popn, params), # not yet implemented
           "SIS" = model_SIS(popn, params),
           "SIR" = model_SIR(popn, params),
           "SEIR" = model_SEIR(popn, params),
           "SIDR" = model_SIDR(popn, params),
           "SEIDR" = model_SIDR(popn, params))
}

