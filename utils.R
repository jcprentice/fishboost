# Format a table as a TSV file for including in the XML file
table_to_tsv_string <- function(x) {
    x |>
        write.table(row.names = FALSE, col.names = FALSE,
                    quote = FALSE, sep = "\t", eol = "\n") |>
        capture.output() |>
        paste0(collapse = "\n") |>
        {\(x) paste0("\n", x, "\n")}()
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
    matrix(x[c(1, 2, 3,
               2, 4, 5,
               3, 5, 6)],
           nrow = 3, ncol = 3,
           dimnames = list(traits, traits))
}


# clamp values
clamp <- function(x, xmin = -Inf, xmax = Inf) {
    max(min(x, xmax), xmin)
}

# message with parameters
msg_pars <- function(x) {
    x[!startsWith(parameter, "Group effect"),
      .(pars = rename_pars(parameter), mean)] |>
        as.data.frame() |>
        capture.output() |>
        paste0(collapse = "\n") |>
        message()
}
