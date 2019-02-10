#---------------------------------------------------#
# Correlation versus importance plot                #
# Bruna Wundervald                                  #        
# January, 2019                                     #        
#---------------------------------------------------#
library(tidyverse)
library(RRF)
library(patchwork)

#setwd("/users/bruna/Project 1")
models <- readRDS("models.rds")
model_paper <- models[[2]] 

# Real gene data ------------------------------------------------------
da <-  data.table::fread("data/PV_Data1_ANON.csv") %>% 
  select(-idPhen) 
dim(da)
nam <- names(da)[-1]

correlations <- readRDS("correlations.rds")

formatted_res <- model_paper %>% 
  purrr::map(RRF::importance) %>%
  setNames(paste0("model", 1:100)) %>% 
  bind_cols() %>% 
  ungroup() %>% 
  mutate(mean = rowSums(.)/100,
         vars = rownames(RRF::importance(model_paper[[1]]))) %>% 
  select(mean, vars) %>% 
  mutate(corr = abs(correlations[abs(correlations) > 0.16]))

formatted_res %>% 
  ggplot(aes(mean, corr)) +
  geom_point(colour = "#e68c7c", size = 2.7, alpha = 0.7) +
  geom_hline(yintercept = 0.5, linetype = "dotted",
             alpha = 0.7, colour = "black", 
             size = 1.3) +
  theme_bw() +
  coord_flip() +
  labs(y = "Absolute values of marginal correlations to y",
       x = 
         expression("Mean importance of the variables"))

