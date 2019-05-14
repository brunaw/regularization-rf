#---------------------------------------------------#
# Regularized Random Forests                        #        
# Bruna Wundervald                                  #        
# April, 2019                                       #        
#---------------------------------------------------#
library(tidyverse) # essential tools
library(RRF) # regularized random forest
library(patchwork)

rmse <- function(model){
  res <- test$y - predict(model, test)
  sqrt(mean(res^2))
}


pars <- function(model){
  length(model$feaSet) 
}


da <- read.table("data/fried_added.txt") 

set.seed(2019)

split <- da %>% 
  mutate(set = ifelse(
    runif(nrow(.)) > 0.75, "test", "train"))

# Train and test set split 
train <- split %>% filter(set == "train") %>% select(-set)
test <- split %>% filter(set == "test") %>% select(-set)

# Using only the correlation as weighting 
corr_fc <- function(var){
  cor(da$y, y = da %>% pull(var), method = "spearman")
}

nam <- names(da)[-6]

corr_vars <- nam %>% 
  purrr::map_dbl(corr_fc)
sort(corr_vars)

selected_vars <- nam[abs(corr_vars) > 0.01] %>% 
  paste(collapse = ' + ')

corr_weight <- corr_vars[abs(corr_vars) > 0.01]

form <- as.formula(paste("y ~ ", selected_vars))

grrf_corr <- RRF(form, data = train,  
                 coefReg =  abs(corr_weight))

pred <- predict(grrf_corr, newdata = test)
mean((pred - test$y)^2) # 4.63

# Selected vars
sel_new <- grrf_corr$importance %>% 
  as.data.frame() %>% 
  mutate(var = rownames(grrf_corr$importance)) %>% 
  arrange(desc(IncNodePurity)) %>% 
  slice(1:10) %>% 
  pull(var) %>% 
  paste(collapse = ' + ')

rrf_new <- randomForest::randomForest(as.formula(
  paste("y ~", sel_new)), data = train)
rrf_new

pred <- predict(rrf_new, newdata = test)
mean((pred - test$y)^2) # 6.024

#-----------------------------------------------------------------------
# Summary of the results NEED TO UPDATE
#-----------------------------------------------------------------------
# Using correlation as the weight
# Vars. selected = 32
# % var explained: 91.24
# MSE = 2.22
# MSE in test set = 4.63
# MSE with the top 10 variables in the test set: 4.005
#-----------------------------------------------------------------------

# First model: sing only the correlation, repeated 250 times
fc_corr_model <- function(i){
  set.seed(i)
  RRF(form, data = train, coefReg =  abs(corr_weight), mtry = 100)
}

system.time(
  models_correlation <- 1:100 %>% 
    purrr::map(fc_corr_model)
)

# Second model: The model from the paper
fc_norm_model <- function(i){
  set.seed(i)
  m1 <- RRF::RRF(form, train, coefReg = 0.8, mtry = 100)
  
  impRF <- m1$importance
  impRF <- impRF[,"IncNodePurity"]    
  imp <- impRF/(max(impRF))
  
  gamma = 0.9
  coefReg <- (1 - gamma)*1 + (gamma*imp )
  set.seed(i)
  RRF(form, data = train, coefReg = coefReg, mtry = 100)
}

system.time(
  model_paper <- 1:100 %>% 
    purrr::map(fc_norm_model)
)

# Third model: using the marginal abs(correlations) *  the usual 
# lambda and gamma as weighting

fc_normc_model <- function(i){
  set.seed(i)
  m1 <- RRF::RRF(form, train, coefReg = 0.9, mtry = 100)
  
  impRF <- m1$importance
  impRF <- impRF[,"IncNodePurity"]    
  
  gamma = 0.9
  n <- data.frame(imp = c(impRF)/(max(c(impRF))),
                  corr = abs(corr_weight)) %>% 
    dplyr::mutate(new_weigth = ifelse(corr > 0.5, 
                                      gamma*imp*corr,
                                      gamma*imp*0.3))
  
  
  set.seed(i)
  RRF(form, data = train, coefReg = n$new_weigth, mtry = 100)
}

system.time(
  model_corre_imp <- 1:100 %>% 
    purrr::map(fc_normc_model)
)

results <- data.frame(
  corr = models_correlation %>% purrr::map_dbl(rmse),
  norm = model_paper %>% purrr::map_dbl(rmse),
  normc = model_corre_imp %>% purrr::map_dbl(rmse))

res_formatted <- results %>% 
  gather("model", "value") %>% 
  mutate(model = forcats::fct_recode(
    model,
    "Using lambda and gamma" = "norm",
    "Correlation" = "corr",
    "Using lambda and gamma * correlation" = "normc"
  ))

head(res_formatted)
p1 <- res_formatted %>% 
  ggplot(aes(x = value, group = model)) +
  geom_density(aes(fill = model), size = 0.7, 
               colour = "grey35", alpha = 0.7) +
  scale_fill_manual(
    values = c("#c03728", "#919c4c", "#f5c04a")) +
  guides(fill = FALSE) +
  #xlim(0.0475, 0.0575)  + 
  labs(x = "Mean squared errors", y = "Density",
       fill = 'Variable weighting', 
       title = "Comparing different weightings for the GRRF",
       subtitle = "Results obtained with 100 runs of each model") +
  theme_bw() 


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
         "Friedman data")  + 
  theme_bw()  +
  theme(legend.position = "bottom",
        legend.direction="vertical")


p1 + p2 + plot_layout(ncol = 1)
