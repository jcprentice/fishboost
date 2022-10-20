# Format a data table as a TSV file for including in the XML file
data_table_to_tsv_string <- function(dt) {
    paste0("\n",
           paste0(capture.output(
               write.table(
                   dt, "", sep = "\t", eol = "\n",
                   quote = FALSE, row.names = FALSE, col.names = FALSE
               )
           ), collapse = "\n"),
           "\n")
}


# Need this safe version of sample() as sometimes length(x) = 1, which leads to
# the "convenience feature", which isn't actually all that convenient.
safe_sample <- function(x, ...) {
    if (length(x) > 1) {
        sample(x, ...)
    } else if (length(x) <= 1) {
        x
    }
}


# slightly nicer way to print an object's size
obsize <- function(x) {
    print(object.size(x), units = "auto")
}


# turn SIRE's 6 values into a covariance matrix
make_cov_matrix <- function(x) {
    traits <- c("sus", "inf", "rec")
    matrix(
        x[c(1, 2, 3,
            2, 4, 5,
            3, 5, 6)],
        nrow = 3,
        ncol = 3,
        dimnames = list(traits, traits)
    )
}


# clamp values
clamp <- function(x, xmin = -Inf, xmax = Inf) {
    max(min(x, xmax), xmin)
}
