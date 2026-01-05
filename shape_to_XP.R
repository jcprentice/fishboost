{
    library(stringr)
    library(purrr)
    library(data.table)
}

shape_to_XP <- function() {
   files <- list.files("datasets", "rds", recursive = TRUE, full.names = TRUE) |>
       str_subset("results") |>
       str_sort(numeric = TRUE)
   
   no_renames <- c()
   
   walk(files, \(f) {
       x <- readRDS(f)
       
       erg <- c("eta", "rho", "gamma")
       rate_pars <- str_c("r_", erg, "_rate")
       
       pars <- expand.grid("r", erg, c("shape", "scale")) |>
           apply(1, str_flatten, "_")
       
       dist_pars <- str_c(c("l", "d", "r"), "p_dist")
       
       write_file <- any(c(rate_pars, pars, dist_pars) %in% names(x$params))
       
       if (write_file) {
           message(str_glue("fixing '{f}'"))
       } else {
           no_renames <<- c(no_renames, f)
           return()
       }
       
       # Rate -> Scale
       idxs <- which(names(x$params) %in% rate_pars)
       
       x$params[idxs] <- map(x$params[idxs], ~ 1 / .x)
       
       names(x$params[idxs]) <- names(x$params[idxs]) |>
           str_replace("rate", "scale")
       
       
       # r_eta -> LP etc.
       idxs <- which(names(x$params) %in% pars)
       
       names(x$params)[idxs] <- names(x$params)[idxs] |>
           str_replace_all(c("r_" = "",
                             "eta" = "LP", "rho" = "DP", "gamma" = "RP"))
       
       # lp_dist -> LP_dist
       idxs <- which(names(x$params) %in% dist_pars)
       
       names(x$params)[idxs] <- names(x$params)[idxs] |>
           str_replace_all(c("lp" = "LP", "dp" = "DP", "rp" = "RP"))
       
       if (write_file) saveRDS(x, f)
   })
   
   if (length(no_renames) > 0) {
       saveRDS(no_renames, "no_renames.rds")
   }
}

shape_to_XP()
