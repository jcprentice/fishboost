get_R0 <- function(popn) {
    has_state <- "state" %in% names(popn)

    if (!has_state) {
        popn[, state := 1L]
    }

    x <- popn[sdp == "progeny",
              .(state, group, gen = map_int(generation, first))]
    x[is.na(gen), gen := 0]

    R0       <- x[, .(R0 = sum(gen == 2) / sum(gen == 1))]$R0
    state_R0 <- x[, .(R0 = sum(gen == 2) / sum(gen == 1)), state]
    group_R0 <- x[, .(R0 = sum(gen == 2) / sum(gen == 1)), .(state, group)]


    message("R0 estimate: ", signif(R0, 3))

    if (!has_state && "state" %in% names(popn)) {
        popn[, state := NULL]
    }

    mget(c("R0", "state_R0", "group_R0"))
}


source_is_loaded <- union(source_is_loaded, "get_R0.R")

