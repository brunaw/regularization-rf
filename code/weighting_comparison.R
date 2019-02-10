#---------------------------------------------------#
# Guided Regularized Random Forests                 #
# Testing differents weights                        #
# Bruna Wundervald                                  #        
# January, 2019                                     #        
#---------------------------------------------------#

library(tidyverse)
library(RRF)
library(patchwork)

#setwd("/users/bruna/Project 1")

# Real gene data ------------------------------------------------------
da <-  data.table::fread("data/PV_Data1_ANON.csv") %>% 
  select(-idPhen) 

dim(da)

# Marginal correlations to the response
corr_fc <- function(var){
  cor(train$log_brd5, y = train %>% pull(var))
}

# Function to calculate the mean squared errors in the test set
msr <- function(model){
  res <- test$log_brd5 - predict(model, test)
  mean(res^2) 
}

# Splitting in train and test data ------------------------------------
set.seed(2019)
split <- da %>% 
  mutate(set = ifelse(
    runif(nrow(.)) > 0.75, "train", "test"))

split %>% 
  janitor::tabyl(set) %>% 
  mutate(percent = scales::percent(percent))

train <- split %>% filter(set == "train") %>% select(-set)
test <- split %>% filter(set == "test") %>% select(-set)

# Applying marginal correlation as the variable weights ----------------
nam <- names(da)[-1]

correlations <- readRDS("correlations.rds")

length(nam[abs(correlations) > 0.16]) # 2626

selected_vars <- nam[abs(correlations) > 0.16] %>% 
  paste(collapse = ' + ')

corr_weight <- correlations[abs(correlations) > 0.16]

form <- as.formula(paste("log_brd5 ~ ", selected_vars))


# First model: using the marginal abs(correlations) as weights

fc_corr_model <- function(i){
  set.seed(i)
  RRF(form, data = train, coefReg =  abs(corr_weight))
}

models_correlation <- 1:100 %>% 
  purrr::map(fc_corr_model)
# Time:  - 10h

# Second model: using the usual lambda and gamma as weighting
fc_norm_model <- function(i){
  set.seed(i)
  m1 <- RRF::RRF(form, train, coefReg = 0.8)
  
  impRF <- m1$importance
  impRF <- impRF[,"IncNodePurity"]    
  imp <- impRF/(max(impRF))
  
  gamma = 0.9
  coefReg <- (1 - gamma)*1 + (gamma*imp )
  set.seed(i)
  RRF(form, data = train, coefReg = coefReg)
}

model_paper <- 1:100 %>% 
  purrr::map(fc_norm_model)


# Third model: using the marginal abs(correlations) *  the usual 
# lambda and gamma as weighting

fc_normc_model <- function(i){
  set.seed(i)
  m1 <- RRF::RRF(form, train, coefReg = 0.8)
  
  impRF <- m1$importance
  impRF <- impRF[,"IncNodePurity"]    
  
  gamma = 0.95
  
  n <- data.frame(imp = c(impRF)/(max(c(impRF))),
                  corr = abs(corr_weight)) %>% 
    dplyr::mutate(new_weigth = ifelse(corr > 0.5, 
                                      (1 - gamma)*1 + (gamma*imp*corr),
                                      (1 - gamma)*1 + (gamma*imp*0.2)))
  
  
  set.seed(i)
  RRF(form, data = train, coefReg = n$new_weigth)
}

model_corre_imp <- 1:100 %>% 
  purrr::map(fc_normc_model)

saveRDS(list(models_correlation, model_paper, 
             model_corre_imp), "models.rds")

# Visualising results --------------------------------------------------
results <- data.frame(
  corr = models_correlation %>% purrr::map_dbl(msr),
  norm = model_paper %>% purrr::map_dbl(msr),
  normc = model_corre_imp %>% purrr::map_dbl(msr))

res_formatted <- results %>% 
  gather("model", "value") %>% 
  mutate(model = forcats::fct_recode(
    model,
    "Using lambda and gamma" = "norm",
    "Correlation" = "corr",
    "Using lambda and gamma * correlation" = "normc"
  ))


p1 <- res_formatted %>% 
  ggplot(aes(x = value, group = model)) +
  geom_density(aes(fill = model), size = 0.7, 
               colour = "grey35", alpha = 0.7) +
  scale_fill_manual(
    values = c("#c03728", "#919c4c", "#f5c04a")) +
  guides(fill = FALSE) +
  xlim(0.0475, 0.0575)  + 
  labs(x = "Mean squared errors", y = "Density",
       fill = 'Variable weighting', 
       title = "Comparing different weightings for the GRRF",
       subtitle = "Results obtained with 100 runs of each model") +
  theme_bw() 


pars <- function(model){
  length(model$feaSet) 
}
  

results_pars <- data.frame(
  corr = models_correlation %>% purrr::map_dbl(pars),
  norm = model_paper %>% purrr::map_dbl(pars),
  normc = model_corre_imp %>% purrr::map_dbl(pars))

res_pars_formatted <- results_pars %>% 
  gather("model", "value") %>% 
  mutate(model = forcats::fct_recode(
    model,
    "Using lambda and gamma" = "norm",
    "Correlation as reg. parameter" = "corr",
    "Using lambda and gamma * correlation" = "normc"
  ))


p2 <- res_pars_formatted %>% 
  ggplot(aes(x = value, group = model)) +
  geom_histogram(aes(fill = model), colour = "grey35",
                 bins = 35,
                 alpha = 0.7, position = "identity") +
  scale_fill_manual(
    guide = guide_legend(),
    values = c("#c03728", "#919c4c", "#f5c04a"),
    labels = c("Correlation as reg. parameter", 
               expression("GRRF:"~lambda~"= 0.8 and"~gamma~"= 0.9"), 
               expression("GRRF: combination of "~lambda~", "~gamma~" and correlation"))) +
  labs(x = "Number of selected variables", y = "Counts",
       fill = 'Variable weighting', 
       caption = 
         "Real gene expression data using the variables with marginal
abs(correlation) to the response > 0.16 (2626)")  + 
  theme_bw()  +
  theme(legend.position = "bottom",
        legend.direction="vertical")

write_rds(list(p1, p2), "plots.rds")

p1 + p2 + plot_layout(ncol = 1)

#-----------------------------------------------------------------------