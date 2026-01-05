{
    library(stringr)
    library(purrr)
    
    source("plotting/plot_covariances.R")
}
    
# dataset <- "fb-final"; scens <- 1:8
dataset <- "fb-test"; scens <- 1:8; scen <- 2;
# dataset <- "fb-final"; scen <- 2L

plts <- map(scens, ~ plot_covariances(dataset, .x))
plt <- plot_covariances(dataset = dataset, scen = scen)

# Slightly scaled for adding to paper
# ggsave(str_glue("{dataset}/gfx/cov-G-{scen}-1.png"),
#        plt$G + theme(text = element_text(size = 10)),
#        width = 24, height = 12, units = "cm")

plt1 <- plot_covariances(dataset = "fb-test", scen = 1L)
plt2 <- plot_covariances(dataset = "fb-test", scen = 8L)

title_plt <- ggplot() + theme_classic()

tplt_ped_cov <- title_plt + labs(title = "Covariances and correlations via pedigree")
tplt_ped_h2  <- title_plt + labs(title = "Heritability and correlations via pedigree")
tplt_grm_cov <- title_plt + labs(title = "Covariances and correlations via GRM")
tplt_grm_h2  <- title_plt + labs(title = "Heritability and correlations via GRM")

plt_ped_cov <- plot_grid(tplt_ped_cov, plt1$G, ncol = 1L, rel_heights = c(0.08, 1))
plt_ped_h2  <- plot_grid(tplt_ped_h2,  plt1$H, ncol = 1L, rel_heights = c(0.08, 1))
plt_grm_cov <- plot_grid(tplt_grm_cov, plt2$G, ncol = 1L, rel_heights = c(0.08, 1))
plt_grm_h2  <- plot_grid(tplt_grm_h2,  plt2$H, ncol = 1L, rel_heights = c(0.08, 1))

gfx_dir <- str_glue("datasets/{dataset}/gfx/cov")

ggsave(str_glue("{gfx_dir}/cov-G-ped.png"),
       plt_ped_cov, width = 9, height = 6)
ggsave(str_glue("{gfx_dir}/cov-h2-ped.png"),
       plt_ped_h2, width = 9, height = 6)
ggsave(str_glue("{gfx_dir}/cov-G-grm.png"),
       plt_grm_cov, width = 9, height = 6)
ggsave(str_glue("{gfx_dir}/cov-h2-grm.png"),
       plt_grm_h2, width = 9, height = 6)

if (FALSE) {
    x <- map(1:10, \(i) {
        plot_covariances(dataset = "fb", scen = 3, itn = i, ci = "hpdi")
    })
}

if (FALSE) {
    x2 <- map(1:8, \(i) {
        plot_covariances(dataset = "sim-test_cov", scen = i, ci = "hpdi")
    })
}
