library(data.table)
library(ggplot2)
library(cowplot)
library(stringr)

dataset <- "fb-parasites4"; scen <- 3

ebvs <- fread(str_glue("data/{dataset}/scen-{scen}-1-out/ebvs.csv"))
fb_data <- readRDS("fb_data/fb_12.rds")

setnames(ebvs, c("id", "inf_g", "inf_e", "sus_g", "sus_e", "tol_g", "tol_e"))
ebvs[, `:=`(trial = as.factor(fb_data$trial),
            donor = as.factor(fb_data$donor),
            sire = fb_data$sire)]

SvIt <- ggplot(ebvs) + geom_point(aes(x = sus_g, y = inf_g, colour = trial))
SvRt <- ggplot(ebvs) + geom_point(aes(x = sus_g, y = tol_g, colour = trial))
IvRt <- ggplot(ebvs) + geom_point(aes(x = inf_g, y = tol_g, colour = trial))

ebvs_trial <- plot_grid(SvIt, SvRt, IvRt)
ggsave("EBVs by trial.pdf", ebvs_trial, width = 10, height = 6)

SvId <- ggplot(ebvs) + geom_point(aes(x = sus_g, y = inf_g, colour = donor))
SvRd <- ggplot(ebvs) + geom_point(aes(x = sus_g, y = tol_g, colour = donor))
IvRd <- ggplot(ebvs) + geom_point(aes(x = inf_g, y = tol_g, colour = donor))

ebvs_donor <- plot_grid(SvId, SvRd, IvRd)
ggsave("EBVs by donor.pdf", ebvs_donor, width = 10, height = 6)

ebvs[, td := str_c("t", trial, "-d", donor)]

SvItd <- ggplot(ebvs) + geom_point(aes(x = sus_g, y = inf_g, colour = td))
SvRtd <- ggplot(ebvs) + geom_point(aes(x = sus_g, y = tol_g, colour = td))
IvRtd <- ggplot(ebvs) + geom_point(aes(x = inf_g, y = tol_g, colour = td))

ebvs_td <- plot_grid(SvItd, SvRtd, IvRtd)
ggsave("EBVs by trial+donor.pdf", ebvs_td, width = 10, height = 6)

