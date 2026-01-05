library(data.table)
library(ggplot2)

x <- readRDS("fb-fes3/results/scen-4-1.rds")
estimated_BVs <- x$estimated_BVs

cols <- c("id", "inf_g", "inf_e", "sus_g", "sus_e", "end_g", "end_e")
setnames(estimated_BVs, cols)
estimated_BVs[, R0 := sus_g + inf_g + end_g]

EBVs1 <- melt(estimated_BVs[, .(id, sus_g, inf_g, end_g, R0)], id.vars = "id")


x <- readRDS("fb-fes3/results/scen-12-1.rds")
estimated_BVs <- x$estimated_BVs

setnames(estimated_BVs, cols)
estimated_BVs[, R0 := sus_g + inf_g + end_g]

EBVs2 <- melt(estimated_BVs[, .(id, sus_g, inf_g, end_g, R0)], id.vars = "id")

EBVs <- rbind(EBVs1, EBVs2, idcol = "set")
EBVs[, set := ifelse(set == 1, "Pedigree", "GRM")]

plt <- ggplot(EBVs, aes(x = value, fill = variable)) +
    geom_density(alpha = 0.5) +
    facet_grid(rows = vars(set),
               cols = vars(variable))

ggsave("gfx/R0_histogram2.png",
       plot = plt, width = 6, height = 6)

