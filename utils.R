# We often just want the first or last letter of a string
str_1st <- function(x) str_sub(x, 1L, 1L)
str_last <- function(x) str_sub(x, -1L)

# split a string into individual characters
str_chars <- function(x) str_split_1(x, "")

# Write sparse matrix so that SIRE can read it
write_sparse_matrix <- function(M, file) {
    dt <- as.data.table(summary(M))
    # Matrix is symmetric, but SIRE needs both halves specified
    dt <- rbind(dt, dt[i != j, .(i = j, j = i, x)])
    setorder(dt, i, j)
    fwrite(dt, file = file,
           sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    dt
}


# Write sparse matrix so that SIRE can read it
write_dt_as_sparse <- function(dt, file) {
    dt <- rbind(dt, dt[i != j, .(i = j, j = i, x)])
    setorder(dt, i, j)
    fwrite(dt, file = file,
           sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    dt
}


# Format a table as a TSV file for including in the XML file
table_to_tsv_string <- function(x, row_names = FALSE, col_names = FALSE,
                                quote = FALSE) {
    x |>
        write.table(row.names = row_names, col.names = col_names,
                    quote = quote, sep = "\t", eol = "\n") |>
        capture.output() |>
        str_flatten("\n") |>
        {\(x) str_c("\n", x, "\n")}()
}

# message a list
capture_message <- function(x) {
    x |> print() |> capture.output() |> str_flatten("\n") |> message()
}


# Need this safe version of sample() as sometimes length(x) = 1, which leads to
# the "convenience feature", which so inconvenient that they even warn you about
# it in the documentation
safe_sample <- function(x, ...) {
    if (length(x) > 1L) {
        sample(x, ...)
    } else if (length(x) <= 1L) {
        x
    }
}

# Merge lists A and B, being careful to not overwrite anything in A from B
safe_merge <- function(A, B) {
    c(A, B[setdiff(names(B), names(A))])
}


# Slightly nicer way to print an object's size
obsize <- function(x) {
    print(object.size(x), units = "auto")
}

nuniq <- function(x) length(unique(x))

# Get all unique characters in a string
uniq_chars <- function(x) str_chars(x) |> unique()


# Clamp values
clamp <- function(x, xmin = -Inf, xmax = Inf) {
    x |> max(xmin) |> min(xmax)
}

# head_tail
ht <- function(x, n = 3L) c(head(x, n), tail(x, n))

# Test if in range
in_range <- function(x, a, b, inc = "ab") {
    switch(inc,
           "ab" = a <= x & x <= b,
           "a" = a <= x & x < b,
           "b" = a < x & x <= b,
           a < x & x < b)
}

# Message with parameters
msg_pars <- function(x) {
    x[!str_starts(parameter, "Group effect"),
      .(pars = rename_pars(parameter),
        mean = signif(mean, 3L),
        hdi95_min = signif(hdi95min, 2L),
        hdi95_max = signif(hdi95max, 2L),
        ESS,
        GR = round(GR, 2))] |>
        as.data.frame() |>
        capture.output() |>
        str_flatten("\n") |>
        message()
}

# Recentering weights
recentre <- function(x, digits = -1L) {
    y <- x - mean(x, na.rm = TRUE)
    if (digits >= 0L) {
        y <- round(y, digits)
    }
    y
}

log_recentre <- function(x, digits = -1L) {
    lx <- log(x)
    y <- lx - mean(lx, na.rm = TRUE)
    if (digits >= 0L) {
        y <- round(y, digits)
    }
    y
}

# Get part of description (used in param_generators)
get_part <- function(x, y) {
    x |>
        str_squish() |>
        str_split_1(", ") |>
        str_subset(y) |>
        str_split_i(" ", 2)
}
