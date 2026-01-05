
reorder_protocol <- function(protocol) {
    # This is necessary to ensure that `expand_priors` comes before any priors,
    # since the priors need to override this.
    if ("expand_priors" %in% names(protocol) && any(str_detect(names(protocol), "prior__"))) {
        protocol <- protocol[append(str_subset(names(protocol), "prior__", negate = TRUE),
                                    str_subset(names(protocol), "prior__"),
                                    after = str_which(names(protocol), "expand_priors"))]
    }
    protocol
}
