set_traits <- function(popn, params) {
    switch(str_to_lower(params$traits_source),
           "grm"       = make_traits_from_grm(popn, params),
           "pedigree"  = make_traits_from_pedigree(popn, params),
           "posterior" = patch_in_traits(popn, params),
           popn)
}

