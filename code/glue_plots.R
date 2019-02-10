library(tidyverse)
library(patchwork)
plots <- readRDS("plots.rds")

p1 <- plots[[1]]
p2 <- plots[[2]]
p1 + p2 + patchwork::plot_layout(ncol = 1)

