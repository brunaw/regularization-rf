#---------------------------------------------------#
# Guided Regularized Random Forests                 #
# Using correlation as variable weighting           #
# Bruna Wundervald                                  #        
# January, 2019                                     #        
#---------------------------------------------------#

library(tidyverse)
library(RRF)

#setwd("/users/bruna/Project 1")

# Real gene data ------------------------------------------------------
da <-  data.table::fread("data/PV_Data1_ANON.csv") %>% 
  select(-idPhen) 

dim(da)

# Marginal correlations to the response
corr_fc <- function(var){
  cor(da$log_brd5, y = train %>% pull(var))
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

corr_vars <- nam[1:20000] %>% 
  purrr::map_dbl(corr_fc)

corr_vars_2 <- nam[20001:length(nam)] %>% 
  purrr::map_dbl(corr_fc)

corr_vars_tot <- c(corr_vars, corr_vars_2)
# Total time: 35 minutes

length(nam[abs(corr_vars_tot) > 0.05] ) # 16021

selected_vars <- nam[abs(corr_vars_tot) > 0.05] %>% 
  paste(collapse = ' + ')

corr_weight <- corr_vars_tot[abs(corr_vars_tot) > 0.05]

form <- as.formula(paste("log_brd5 ~ ", selected_vars))

grrf_corr <- RRF(form, data = train,  
                 coefReg =  abs(corr_weight))

# Modeling time: 11 minutes
sum(grrf_corr$importance > 0) # 96 important variables 
pred <- predict(grrf_corr, newdata = test)
mean((pred - test$log_brd5)^2) 

# Selected vars
sel_new <- grrf_corr$importance %>% 
  as.data.frame() %>% 
  mutate(var = rownames(grrf_corr$importance)) %>% 
  arrange(desc(IncNodePurity)) %>% 
  slice(1:10) %>% 
  pull(var) %>% 
  paste(collapse = ' + ')

rrf_new <- randomForest::randomForest(as.formula(
  paste("log_brd5 ~", sel_new)), data = train)
rrf_new

pred <- predict(rrf_new, newdata = test)
mean((pred - test$log_brd5)^2) 

#-----------------------------------------------------------------------
# Summary of the results
#-----------------------------------------------------------------------
# Using correlation as the weight
# Vars. selected = 100
# % var explained: 43.4
# MSE = 0.051
# MSE in test set = 0.050
# MSE with the top 10 variables in the test set: 0.054
#-----------------------------------------------------------------------
# Pre-selecting variables using the correlation criteria
# Vars. selected = 70
# % var explained: 59.63%
# MSE = 0.037
# MSE in test set = 0.0531
# MSE with the top 10 variables in the test set: 0.057
#-----------------------------------------------------------------------
# Using correlation * regularizers as the weight  
# Vars. selected = 71
# % var explained: 46.7
# MSE = 0.048
# MSE in test set = 0.0538
# MSE with the top 10 variables in the test set: 0.057
#-----------------------------------------------------------------------

# Usual method
m1 <- RRF::RRF(form, train, coefReg = 0.8)

impRF <- m1$importance
impRF <- impRF[,"IncNodePurity"]              # get the importance score
imp <- impRF/(max(impRF))

gamma = 0.9
coefReg <- (1 - gamma)*1 + (gamma*imp )
grrf <- RRF(form, data = train, coefReg = coefReg)

sum(grrf$importance > 0) # 70 important variables 

pred <- predict(grrf, newdata = test)
mean((pred - test$log_brd5)^2) 

# Selected vars
sel_new <- grrf$importance %>% 
  as.data.frame() %>% 
  mutate(var = rownames(grrf$importance)) %>% 
  arrange(desc(IncNodePurity)) %>% 
  slice(1:10) %>% 
  pull(var) %>% 
  paste(collapse = ' + ')

rrf_new <- randomForest::randomForest(as.formula(
  paste("log_brd5 ~", sel_new)), data = train)
rrf_new

pred <- predict(rrf_new, newdata = test)
mean((pred - test$log_brd5)^2) 


# Correlation * usual criteria 
grrf_and_corr <- RRF(form, data = train, 
                     coefReg = coefReg * abs(corr_weight))

sum(grrf_and_corr$importance > 0) # 71 important variables 

pred <- predict(grrf_and_corr, newdata = test)
mean((pred - test$log_brd5)^2) 

# Selected vars
sel_new <- grrf$importance %>% 
  as.data.frame() %>% 
  mutate(var = rownames(grrf$importance)) %>% 
  arrange(desc(IncNodePurity)) %>% 
  slice(1:10) %>% 
  pull(var) %>% 
  paste(collapse = ' + ')

rrf_corr_new <- randomForest::randomForest(as.formula(
  paste("log_brd5 ~", sel_new)), data = train)
rrf_corr_new

pred <- predict(rrf_corr_new, newdata = test)
mean((pred - test$log_brd5)^2) 
