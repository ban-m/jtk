library(tidyverse)
loadNamespace("cowplot")

filename <- "./result/obcx/jtk_result/clustering_bench/two_cluster.tsv" # nolint
col_names <- c("key", "seed", "length", "time",
              "rand", "ARI", "coverage", "error") # nolint
dataset <- read_tsv(filename, col_names = col_names)

g <- dataset %>%
    filter(coverage %in% c(10, 20, 30)) %>%
    filter(error != 0.05) %>%
    filter(length %in% c(1000, 2000, 4000)) %>%
    ggplot() +
    geom_boxplot(aes(x = factor(coverage), y = rand)) +
    geom_hline(aes(yintercept = 1 - error), color = "red") +
    facet_grid(error ~ length) +
    labs(x = "Coverage (per haploid)", y = "Rand Index") +
    cowplot::theme_cowplot()

cowplot::ggsave2("./result/plots/two_cluster.png", plot = g,
    width = 15, height = 10, units = "cm", dpi = 500)

filename <- "./result/obcx/jtk_result/clustering_bench/four_cluster.tsv" # nolint
col_names <- c("key", "seed", "length", "time", "rand", "ARI", "coverage", "error") # nolint
dataset <- read_tsv(filename, col_names = col_names)

g <- dataset %>%
    filter(length %in% c(1000, 2000, 4000)) %>%
    filter(error != 0.05) %>%
    filter(coverage %in% c(10, 20, 30)) %>%
    ggplot() +
    geom_boxplot(aes(x = factor(coverage), y = rand)) +
    geom_hline(aes(yintercept = 0.8), color = "red") +
    facet_grid(error ~ length) +
    labs(x = "Coverage (per haploid)", y = "Rand Index") +
    cowplot::theme_cowplot()


cowplot::ggsave2("./result/plots/four_cluster.png", plot = g,
    width = 15, height = 10, units = "cm", dpi = 500)
