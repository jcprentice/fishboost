# Adds a check to avoid re-sourcing files all over the place, also lets you keep
# track of everything that's been loaded already.

safe_source <- function(...) {
    if ("src_loaded" %notin% names(.GlobalEnv)) {
        .GlobalEnv$src_loaded <- c()
    }

    for (s in list(...)) {
        if (s %notin% src_loaded) {
            message("- sourcing '", s, "'")
            source(s)
            .GlobalEnv$src_loaded <- union(.GlobalEnv$src_loaded, s)
        }
    }
}
