# Adds a check to avoid re-sourcing files all over the place, also lets you keep
# track of everything that's been loaded already.

safe_source <- function(...) {
    if ("source_is_loaded" %notin% names(.GlobalEnv)) {
        .GlobalEnv$source_is_loaded <- c()
    }

    for (s in list(...)) {
        if (s %notin% source_is_loaded) {
            message("- sourcing '", s, "'")
            source(s)
            .GlobalEnv$source_is_loaded <- union(.GlobalEnv$source_is_loaded, s)
        }
    }
}
