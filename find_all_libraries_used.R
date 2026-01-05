if (!require(tidyverse)) install.packages("tidyverse")

library(tidyverse)

libs <- list.files(pattern = "\\.R$",
                   recursive = TRUE,
                   full.names = TRUE) |>
    str_subset("find_all_libraries_used",
               negate = TRUE) |>
    map(~ readLines(.x) |>
            str_subset("library") |>
            str_remove_all(" |library\\(|\\)|#[^\n]*")) |>
    flatten_chr() |>
    unique() |>
    str_subset(".+") |>
    str_sort()

install.packages(libs)
