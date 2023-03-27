library(tidyverse)
loadNamespace("cowplot")

filename <- "./result/obcx/jtk_result/clustering_bench/two_cluster.tsv" # nolint
col_names <- c(
    "key", "seed", "length", "time",
    "rand", "ARI", "coverage", "error"
) # nolint
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

cowplot::ggsave2("./result/plots/two_cluster.png",
    plot = g,
    width = 15, height = 10, units = "cm", dpi = 500
)

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


cowplot::ggsave2("./result/plots/four_cluster.png",
    plot = g,
    width = 15, height = 10, units = "cm", dpi = 500
)


filename <- "./result/obcx/jtk_result/num_of_contigs.tsv"
dataset <- read_tsv(filename, col_names = c("type", "length", "id", "contigs"))

g <- dataset %>%
    filter(length > 15000) %>%
    ggplot() +
    geom_point(aes(x = length, y = contigs), alpha = 0.2, size = 4) +
    facet_grid(type ~ .) +
    cowplot::theme_cowplot() +
    labs(
        x = "Maximum length of the reads (bp)",
        y = "Number of phased contigs"
    )

cowplot::ggsave2("./result/plots/haploid_contigs.png",
    plot = g,
    width = 15, height = 10, units = "cm", dpi = 500
)
