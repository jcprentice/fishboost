# setup makes multiple changes to params, which either need to be provided as an
# arg to make_parameters(), or carefully patched afterwards.

apply_setup <- function(setup, params) {
    switch(setup,
           "fb_12_drop71" = {
               nsires <- 28L; ndams <- 25L; nprogeny <- 1750L;
               dpsire <- 1L; ppdam <- 72L;
               ngroups <- 70L; I0 <- 5L;
           }, "fb_1_drop71" = {
               nsires <- 14L; ndams <- 14L; nprogeny <- 875L;
               dpsire <- 1L; ppdam <- 72L;
               ngroups <- 35L; I0 <- 5L;
           }, "fb_2_drop71" = {
               nsires <- 17L; ndams <- 14L; nprogeny <- 875L;
               dpsire <- 1L; ppdam <- 72L;
               ngroups <- 35L; I0 <- 5L;
           }, "fb_12" = {
               nsires <- 29L; ndams <- 25L; nprogeny <- 1775L;
               dpsire <- 1L; ppdam <- 72L;
               ngroups <- 71L; I0 <- 5L;
           }, "fb_1" = {
               nsires <- 14L; ndams <- 14L; nprogeny <- 875L;
               dpsire <- 1L; ppdam <- 72L;
               ngroups <- 35L; I0 <- 5L;
           }, "fb_2" = {
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

    vars <- c("setup", "nsires", "ndams", "nprogeny", "nparents", "ntotal",
              "dpsire", "ppdam", "group_size", "ngroups", "I0")

    params[vars] <- mget(vars)

    params
}
