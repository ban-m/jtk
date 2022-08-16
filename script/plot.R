library(tidyverse)
loadNamespace("cowplot")

filename <- "./result/clustering_bench/two_cluster.tsv"
col_names <- c("key", "seed", "time", "rand", "ARI", "coverage", "error")
dataset <- read_tsv(filename, col_names = col_names)

g <- dataset %>%
    filter(error != 0.02) %>%
    ggplot() +
    geom_boxplot(aes(x = factor(coverage), y = rand)) +
    facet_wrap(vars(error)) +
    labs(x = "Coverage (per haploid)", y = "Rand Index") +
    cowplot::theme_cowplot()

cowplot::ggsave2("./result/plots/two_cluster.png", plot = g,
    width = 20, height = 10, units = "cm", dpi = 500)

filename <- "./result/clustering_bench/four_cluster.tsv"
col_names <- c("key", "seed", "time", "rand", "ARI", "coverage", "error")
dataset <- read_tsv(filename, col_names = col_names)

g <- dataset %>%
    filter(error != 0.02) %>%
    ggplot() +
    geom_boxplot(aes(x = factor(coverage), y = rand)) +
    facet_wrap(vars(error)) +
    labs(x = "Coverage (per haploid)", y = "Rand Index") +
    cowplot::theme_cowplot()


cowplot::ggsave2("./result/plots/four_cluster.png", plot = g,
    width = 20, height = 10, units = "cm", dpi = 500)
