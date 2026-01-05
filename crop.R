crop <- function(popn, trait = "r0", prop = 0.2) {
    # trait <- "r0"; prop <- 0.2
    
    # Work with a copy of just the progeny so we don't mess with popn (since
    # it's a DT, any changes inside the function would affect the DT outside).
    
    popn2 <- popn[sdp == "progeny"]
    
    if (trait == "none") {
        return(popn2)
    }
    
    pop_size <- popn2[, .N]
    group_sizes <- popn2[, .N, group]
    
    # Sort by trait, label top p with cull = TRUE, and resort by id
    setorderv(popn2, trait)
    popn2[, cull := .I > (1 - prop) * .N]
    setorder(popn2, id)
    
    # These are the traits whose values we want to copy
    sit <- c("sus", "inf", "tol")
    cols <- intersect(names(popn2),
                      c(str_c(sit, "_g"), str_c(sit, "_e"), sit))
    
    # Sample IDs from the unculled set
    num_culled <- popn2[cull == TRUE, .N]
    popn2b <- popn2[cull == FALSE, ..cols][sample(num_culled, replace = TRUE)]
    
    # And replace the culled IDs with the sampled set
    popn2[cull == TRUE, (cols) := popn2b]
    
    popn2[, cull := NULL]
    
    # Finally, add the parents back in
    popn3 <- rbind(popn[sdp != "progeny"], popn2)
    
    popn3
}
