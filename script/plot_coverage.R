#!/bin/env R
library(tidyverse)
loadNamespace("cowplot")
args <- commandArgs(trailingOnly = TRUE)
print(args)
cov_dataset <- read_tsv(args[1])
plot <- cov_dataset %>% ggplot() +
    geom_line(aes(x = Position, y = Coverage)) + facet_grid(Contig ~ .) +
    labs(x = "Position(bp)", y = "Coverage") +
    cowplot::theme_cowplot()
cowplot::ggsave2(filename = args[2], plot = plot)
dev.off()