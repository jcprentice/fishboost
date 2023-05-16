source("plotting/covariances.R")

data_set <- "fb-parasites4"
scenario <- 3
replicates <- c(1, 3, 4)
plt <- plot_covariances(data_set = data_set, scenario = scenario, ci = "hdpi")

ggsave(glue("gfx/{data_set}/cov_G-{scenario}.png"),
       plt$G + theme(text = element_text(size = 10)),
       width = 24, height = 12, units = "cm")


if (FALSE) {
    x <- list()
    for (i in 1:10) {
        x[[i]] = plot_covariances(data_set = "fb", scenario = 3, ci = "hpdi", replicate = i)
    }
}

if (FALSE) {
    x2 <- list()
    for (i in 1:8) {
        x2[[i]] = plot_covariances(data_set = "sim-test_cov", scenario = i, ci = "hpdi")
    }
}

