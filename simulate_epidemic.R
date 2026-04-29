simulate_epidemic <- function(popn, params) {
    model_type <- params$model_type
    f <- str_glue("model_{model_type}")
    sf <- str_glue("models/{f}.R")
    if (!file.exists(sf)) {
        stop(str_glue("Sorry, {model_type} not yet implemented!"))
    }
    source(sf)
    get(f)(popn, params)
}

