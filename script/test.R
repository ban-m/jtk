library(tidyverse)
loadNamespace("cowplot")

xlab <- "Length (bp) sorted in ascending order"
ylab <- "Accumurated length / total length (between 0 and 1)"
cum_ttl <- "The read length when the y-axis is 0.5 corresponds to the N50."

hg002 <- "./result/obcx/jtk_result/reads_metric/hg002_reads.tsv"
b080 <- "./result/pg/jtk/reads_metric/lengths.tsv"

hg002_datasets <- read_tsv(hg002, col_names = c("quality", "length", "type"))

g <- hg002_datasets %>%
    filter(length < 200000) %>%
    ggplot() +
    geom_histogram(aes(x = length, y = after_stat(density))) +
    facet_wrap(vars(type)) +
    cowplot::theme_cowplot() +
    labs(x = "Read Length (bp)", y = "Density (%)")

cowplot::ggsave2("./result/plots/histogram_hg002.png",
    plot = g,
    width = 25, height = 20, units = "cm", dpi = 250
)


g <- hg002_datasets %>%
    group_by(type) %>%
    group_modify(~ .x %>%
        arrange(length) %>%
        mutate(cumlen = cumsum(length) / sum(length))) %>%
    ungroup() %>%
    filter(length < 200000) %>%
    ggplot() +
    geom_line(aes(x = length, y = cumlen)) +
    geom_hline(yintercept = 0.5, color = "red") +
    facet_wrap(vars(type)) +
    labs(x = xlab, y = ylab, title = cum_ttl) +
    cowplot::theme_cowplot()

cowplot::ggsave2("./result/plots/cumlen_hg002.png",
    plot = g,
    width = 25, height = 20, units = "cm", dpi = 250
)


b080_datasets <- read_tsv(b080, col_names = c("quality", "length", "type"))
g <- b080_datasets %>%
    filter(length < 200000) %>%
    ggplot() +
    geom_histogram(aes(x = length, y = after_stat(density))) +
    facet_wrap(vars(type)) +
    cowplot::theme_cowplot() +
    labs(x = "Read Length (bp)", y = "Density (%)")

cowplot::ggsave2("./result/plots/histogram_b080.png",
    plot = g,
    width = 25, height = 20, units = "cm", dpi = 250
)


g <- b080_datasets %>%
    group_by(type) %>%
    group_modify(~ .x %>%
        arrange(length) %>%
        mutate(cumlen = cumsum(length) / sum(length))) %>%
    ungroup() %>%
    filter(length < 200000) %>%
    ggplot() +
    geom_line(aes(x = length, y = cumlen)) +
    geom_hline(yintercept = 0.5, color = "red") +
    facet_wrap(vars(type)) +
    labs(x = xlab, y = ylab, title = cum_ttl) +
    cowplot::theme_cowplot()

cowplot::ggsave2("./result/plots/cumlen_b080.png",
    plot = g,
    width = 25, height = 20, units = "cm", dpi = 250
)
