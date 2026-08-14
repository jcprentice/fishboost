{
    library(stringr)
    safe_source("utils.R")
}

# Map "sit" -> c("s", "i", "t"), "si_tt" -> c("s", "i", "t", "t")
get_l1 <- function(tr, traits) {
    x <- names(traits[str_chars(tr)])
    x[!is.na(x)]
}

# Map "sit" -> c("sus", "inf", "tol"), "si_tt" -> c("sus", "inf", "tol", "tol")
get_ltrait <- function(tr, traits) {
    x <- unname(traits[str_chars(tr)])
    x[!is.na(x)]
}

# Map "sit" -> c(1, 2, 5), "si_tt" -> c(1, 2, 5, 5)
get_lidx <- function(tr, traits) {
    x <- match(str_chars(tr), names(traits))
    x[!is.na(x)]
}

# Map "sit" -> c("s", "i", "t"), "si_tt" -> c("s", "i", "d", "t")
get_t1 <- function(tr, traits) {
    x <- names(traits[str_chars(tr)])
    names(traits)[!is.na(x)]
}

# Map "sit" -> c("sus", "inf", "tol"), "si_tt" -> c("sus", "inf", "det", "tol")
get_trait <- function(tr, traits) {
    x <- unname(traits[str_chars(tr)])
    unname(traits[!is.na(x)])
}

# Map "sit" -> c(1, 2, 5), "si_tt" -> c(1, 2, 4, 5)
get_idx <- function(tr, traits) {
    t1 <- str_chars(tr)
    if (length(t1) == length(traits)) {
        x <- match(t1, names(traits))
        seq_along(traits)[!is.na(x)]
    } else {
        match(t1, names(traits))
    }
}

if (FALSE) {
    map(test, get_trait, traits) |> setNames(test)

    micrombenchmark::microbenchmark(
        a = get_l1(tr, traits),
        b = get_ltrait(tr, traits),
        c = get_lidx(tr, traits),
        d = get_t1(tr, traits),
        e = get_trait(tr, traits),
        f = get_idx(tr, traits),
        times = 1e5,
        setup = {tr <- sample(test, 1L)}
    )
}


source_is_loaded <- union(source_is_loaded, "map_idxs.R")

