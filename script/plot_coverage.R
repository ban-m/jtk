#!/bin/env R
library(tidyverse)
loadNamespace("cowplot")
args <- commandArgs(trailingOnly = TRUE)
print(args)
cov_dataset <- read_tsv(args[1])
outfile <- args[2]
plot <- cov_dataset %>% ggplot() +
    geom_line(aes(x = position, y = coverage)) + facet_grid(contig ~ .) +
    labs(x = "Position(bp)", y = "Coverage") +
    scale_y_continuous(breaks = seq(0, 60, length = 5), limits = c(0, 60)) +
    cowplot::theme_cowplot()
cowplot::ggsave2(filename = outfile, plot = plot)
dev.off()
