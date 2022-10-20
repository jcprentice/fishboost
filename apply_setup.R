# setup makes multiple changes to params, which either need to be provided as an
# arg to make_parameters(), or carefully patched afterwards.

apply_setup <- function(setup, params) {
    switch(setup,
           "fishboost" = {
               nsires <- 29L; ndams <- 25L; nprogeny <- 1800L;
               dpsire <- 1L; ppdam <- 72L;
               ngroups <- 72L; I0 <- 5L;
           }, "fb1" = {
               nsires <- 14L; ndams <- 14L; nprogeny <- 900L;
               dpsire <- 1L; ppdam <- 72L;
               ngroups <- 36L; I0 <- 5L;
           }, "fb2" = {
               nsires <- 18L; ndams <- 14L; nprogeny <- 900L;
               dpsire <- 1L; ppdam <- 72L;
               ngroups <- 36L; I0 <- 5L;
           }, "chris" = {
               nsires <- 100L; ndams <- 2000L; nprogeny <- 2000L;
               dpsire <- 20L; ppdam <- 1L;
               ngroups <- 200L; I0 <- 1L;
           }, "small" = {
               nsires <- 3L; ndams <- 6L; nprogeny <- 12L;
               dpsire <- 2L; ppdam <- 2L;
               ngroups <- 4L; I0 <- 1L;
           }, {
               nsires <- 3L; ndams <- 6L; nprogeny <- 12L;
               dpsire <- 2L; ppdam <- 2L;
               ngroups <- 4L; I0 <- 1L;
           }
    )

    # derived numbers
    nparents <- nsires + ndams
    ntotal <- nprogeny + nparents
    group_size <- nprogeny / ngroups

    new_vals <- mget(c("setup", "nsires", "ndams", "nprogeny", "nparents",
                       "ntotal", "dpsire", "ppdam", "group_size", "ngroups",
                       "I0"))

    for (n in names(new_vals)) {
        params[[n]] <- new_vals[[n]]
    }

    params
}
