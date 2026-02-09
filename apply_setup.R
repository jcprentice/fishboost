#' Apply `setup` to `params`
#'
#' @descrition Modify `params` according to `setup.` This is only necessary if
#'   the `popn` data.table needs to be constructed from scratch, otherwise the
#'   values set are never used.
#'
#' @param setup A string describing `popn`'s family structure
#' @param params A params list
#'
#' @returns A new modified params list

apply_setup <- function(setup, params) {
    switch(setup,
           "fb_12_rpw" = {
               nsires <- 28L; ndams <- 25L; nprogeny <- 1750L;
               dpsire <- 1L; ppdam <- 72L;
               ngroups <- 70L; I0 <- 5L;
           }, "fb_1_rpw" = {
               nsires <- 14L; ndams <- 14L; nprogeny <- 875L;
               dpsire <- 1L; ppdam <- 72L;
               ngroups <- 35L; I0 <- 5L;
           }, "fb_2_rpw" = {
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
