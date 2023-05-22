library(data.table)
library(ggplot2)
library(cowplot)

ebvs <- fread("data/fb-parasites4/scen-3-1_out/ebvs.csv")
fb_data <- readRDS("fb_data/fb_data12.rds")

setnames(ebvs, c("id", "inf_G", "inf_E", "sus_G", "sus_E", "rec_G", "rec_E"))
ebvs[, `:=`(trial = as.factor(fb_data$trial),
            donor = as.factor(fb_data$donor),
            sire = fb_data$sire)]

SvIt <- ggplot(ebvs) + geom_point(aes(x=sus_G, y=inf_G, colour = trial))
SvRt <- ggplot(ebvs) + geom_point(aes(x=sus_G, y=rec_G, colour = trial))
IvRt <- ggplot(ebvs) + geom_point(aes(x=inf_G, y=rec_G, colour = trial))

ebvs_trial <- plot_grid(SvIt, SvRt, IvRt)
ggsave("EBVs by trial.pdf", ebvs_trial, width = 10, height = 6)

SvId <- ggplot(ebvs) + geom_point(aes(x=sus_G, y=inf_G, colour = donor))
SvRd <- ggplot(ebvs) + geom_point(aes(x=sus_G, y=rec_G, colour = donor))
IvRd <- ggplot(ebvs) + geom_point(aes(x=inf_G, y=rec_G, colour = donor))

ebvs_donor <- plot_grid(SvId, SvRd, IvRd)
ggsave("EBVs by donor.pdf", ebvs_donor, width = 10, height = 6)

ebvs[, td := paste0("t", trial, "-d", donor)]

SvItd <- ggplot(ebvs) + geom_point(aes(x=sus_G, y=inf_G, colour = td))
SvRtd <- ggplot(ebvs) + geom_point(aes(x=sus_G, y=rec_G, colour = td))
IvRtd <- ggplot(ebvs) + geom_point(aes(x=inf_G, y=rec_G, colour = td))

ebvs_td <- plot_grid(SvItd, SvRtd, IvRtd)
ggsave("EBVs by trial+donor.pdf", ebvs_td, width = 10, height = 6)
