library(tidyverse)
dataset <- read_tsv("./result/two_haplotypes.tsv")
g <- dataset %>%
    ggplot() + geom_histogram(mapping=aes(x=Dist1, fill=factor(Answer)))
ggsave("./png/two_haplotype_hist.png", g)

g <- dataset %>%
    ggplot() + geom_point(mapping=aes(x=Dist1, y= Dist2, color=factor(Answer)))
ggsave("./png/two_haplotype_point.png", g)

g <- dataset %>%
    ggplot() + geom_histogram(mapping=aes(x=Dist1-Dist2, fill=factor(Answer)))
ggsave("./png/two_haplotype_diff_hist.png", g)
