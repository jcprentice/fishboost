{
    library(purrr)
    library(data.table)
    library(stringr)
}

rep_fb <- function(n = 2, base = "fb_12") {
    # n <- 3; base <- "fb_12_drop71"
    
    x <- str_glue("fb_data/{base}.rds") |>
        readRDS() |>
        split(by = "trial") |>
        setNames(c("parents", "t1", "t2"))
    
    ngroups <- x[-1] |>
        map(~ .x[, length(unique(group))]) |>
        unname() |>
        rep(each = n) |>
        cumsum() |>
        append(0, after = 0)
    
    xs <- c(list(x$parents),
            rep(list(x$t1), n),
            rep(list(x$t2), n)) |>
        rbindlist(idcol = "idcol")
    
    walk(seq(2, 2 * n + 1), \(i) {
        xs[idcol == i, group := group - min(group) + ngroups[[i - 1]] + 1]
    })
    
    xs[, `:=`(id = .I, idcol = NULL)]
    
    baseN <- str_replace(base, "fb", str_glue("fb{n}"))
    saveRDS(xs, str_glue("fb_data/{baseN}.rds"))
}

rep_fb(2)
