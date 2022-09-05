library(tidyverse)
loadNamespace("cowplot")

filename <- "lap.tsv"
dataset <- read_tsv(filename, col_names = FALSE)
dataset <- as.matrix(dataset)
heatmap(x = dataset, Rowv = NA, symm = TRUE, scale = "none")


filename <- "341.tsv" # nolint

dataset <- read_tsv(filename, col_names = FALSE)
len <- dataset %>% pull(X1) %>% max()
idmap <- setNames(object = seq(from = 0, to = len),
 nm = dataset %>% filter(X2 == -2) %>% arrange(X3) %>% pull(X1))
dataset$newID <- idmap[as.character(dataset$X1)]

g <- dataset %>%
    ggplot() +
    geom_tile(aes(x = X2, y = newID, fill = X3)) +
    scale_fill_gradient2(low = "blue", high = "orange")

temp <- dataset %>%
filter(0 <= X2) %>%
    select(-newID) %>%
    mutate(X2 = paste0("Pos", X2)) %>%
    pivot_wider(names_from = X2, values_from = X3)

temp %>%
    summarize(across(starts_with("Pos"), ~ sum(.x[0 < .x])))
temp %>%
    summarize(across(starts_with("Pos"), ~ sum(0 < .x)))


matrix_data <- dataset %>% select(-X1, -X2) %>% as.matrix()


col_names <- c("ID", "Contig", "Length", "Position")
correct_aln <- read_tsv("answer.tsv", col_names = col_names)
proposed_aln <- read_tsv("prop.tsv", col_names = c("ID", "assembly"))

joined_data <- full_join(correct_aln, proposed_aln, by = "ID")


diff_data <- read_tsv("pj73.diff")
g <- diff_data %>%
    ggplot() + geom_histogram(aes(x = Qpos), bins = 60) +
    facet_grid(Type ~ Query) +
    cowplot::theme_cowplot() +
    labs(x = "Position in the assembly (bp)", y = "# of in/del/mism")

cowplot::ggsave2("./result/plots/pj73_diff.png", g)


g <- diff_data %>%
    filter(Size > 3) %>%
    ggplot() + geom_histogram(aes(x = Qpos), bins = 60) +
    facet_grid(Type ~ Query) +
    cowplot::theme_cowplot() +
    labs(x = "Position in the assembly (bp)", y = "# of in/del/mism")

cowplot::ggsave2("./result/plots/pj73_diff_above3.png", g)
