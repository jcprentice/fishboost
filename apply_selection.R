library(data.table)

#' Replicate selection by removing and resampling
#'
#' @description `apply_selection()` takes a population data.table and a
#'   parameters list, removes the top x% and replaces them with individuals
#'   sampled from the rest of the population
#'
#' @param popn A population data.table
#' @param params A list of parameters
#'
#' @returns A new population file

apply_selection <- function(popn, params) {
    select_on  <- params$select_on %||% "none"
    select_top <- params$select_top %||% 0.2
    traits_g   <- str_c(params$model_traits, "_g")
    traits_e   <- str_c(params$model_traits, "_e")

    if (is.null(select_on) || select_on == "none") {
        return(popn)
    }

    popn2 <- copy(popn)

    # We don't have an r0_g column, so temporarily create one
    if (select_on == "r0") {
        popn2[, r0_g := sus_g + inf_g + tol_g]
    }
    ids <- popn2[sdp == "progeny"] |>
        _[order(- get(str_c(select_on, "_g")))] |>
        _[seq(select_top * .N), id] |>
        sort()
    if (select_on == "r0") {
        popn2[, r0_g := NULL]
    }

    sample_ids <- setdiff(popn2[sdp == "progeny", id], ids)

    cols <- intersect(c("sire", "dam", traits_g, traits_e), names(popn))
    tmp <- popn2[sample(sample_ids, length(ids), replace = TRUE), ..cols]
    popn2[ids, (cols) := tmp]

    popn2
}
