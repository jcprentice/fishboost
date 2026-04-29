simulate_epidemic <- function(popn, params) {
    model_type <- params$model_type
    f <- str_glue("model_{model_type}")
    source(str_glue("models/{f}.R"))
    get(f)(popn, params)
}

