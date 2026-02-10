library(stringr)
library(purrr)

# For BICI config files, modify them to increment the seed by 100 if they fail
# to run properly
inc_seed <- function(params) {
    f <- params$config

    lines <- readLines(f)
    n1 <- str_which(lines, "inf.*seed")
    line <- lines[[n1]]
    parts <- str_split_1(line, " ")
    n2 <- str_which(parts, "seed")
    part <- parts[[n2]]
    seed <- str_split_i(part, "=", 2) |> as.numeric()
    parts[[n2]] <- str_c("seed=", seed + 100)
    lines[[n1]] <- str_flatten(parts, " ")
    
    writeLines(lines, f)
}
