library(tidyverse)
loadNamespace("cowplot")

filename <- "40.tsv" # nolint

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
