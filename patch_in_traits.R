patch_in_traits <- function(pedigree, params) {
    {
        data_set <- params$patch_data_set
        scenario <- params$patch_scenario
        traitnames <- params$traitnames
    }
    
    f <- glue("data/{data_set}/scen-{scenario}-1_out/ebvs.csv")
    ebvs <- fread(f)
    
    if (nrow(pedigree) != nrow(ebvs)) {
        message(glue("Inconsistent sizes for EBVs ({x}) and pedigree ({y}), generating new data",
                     x = nrow(pedigree), y = nrow(ebvs)))
        traits <- make_traits_from_pedigree(pedigree, params)
        return(traits)
    }
    
    message(glue("Making traits by copying EBVs from {data_set} / scenario {scenario} ..."))
    
    new_names <- copy(names(ebvs))
    new_names <- sub("_a EBV", "_BV", new_names)
    new_names <- sub("_e EBV", "_EV", new_names)
    setnames(ebvs, new_names)
    
    if (!"s_BV" %in% new_names) ebvs[, `:=`(s_BV = 0, s_EV = 0)]
    if (!"i_BV" %in% new_names) ebvs[, `:=`(i_BV = 0, i_EV = 0)]
    if (!"r_BV" %in% new_names) ebvs[, `:=`(r_BV = 0, r_EV = 0)]
    
    setcolorder(ebvs, c("id", "s_BV", "s_EV", "i_BV", "i_EV", "r_BV", "r_EV"))
    
    ebvs[, `:=`(s_ = exp(s_BV + s_EV),
                i_ = exp(i_BV + i_EV),
                r_ = exp(r_BV + r_EV),
                latency = 1,
                detectability = 1)]
    
    new_names <- copy(names(ebvs))
    new_names <- sub("s_", "susceptibility_", new_names)
    new_names <- sub("i_", "infectivity_", new_names)
    new_names <- sub("r_", "recoverability_", new_names)
    new_names <- sub("_$", "", new_names)
    setnames(ebvs, new_names)
    
    traits <- merge(pedigree, ebvs)
    
    traits
}