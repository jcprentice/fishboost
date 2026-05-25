get_R0 <- function(popn) {
    has_state <- "state" %in% names(popn)

    if (!has_state) {
        popn[, state := 1]
    }

    R0 <- popn[sdp == "progeny",
               .(R0 = sum(generation == 2, na.rm = TRUE) /
                     sum(generation == 1, na.rm = TRUE))]$R0
    state_R0 <- popn[sdp == "progeny",
                     .(R0 = sum(generation == 2, na.rm = TRUE) /
                           sum(generation == 1, na.rm = TRUE)),
                     state]
    group_R0 <- popn[sdp == "progeny",
                     .(R0 = sum(generation == 2, na.rm = TRUE) /
                           sum(generation == 1, na.rm = TRUE)),
                     .(state, group)]



    message("R0 estimate: ", signif(R0, 3))

    if (!has_state) {
        popn[, state := NULL]
    }

    mget(c("R0", "state_R0", "group_R0"))
}

