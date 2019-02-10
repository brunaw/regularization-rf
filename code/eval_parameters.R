#--------------------------------------------------------#
# Regularized Random Forests                             #    
# Code for the evaluation of regularization parameters   #
# Bruna Wundervald                                       #        
# January, 2019                                          #        
#--------------------------------------------------------#
library(tidyverse)
library(RRF)
library(patchwork)

lambdas <- seq(0.01, 0.99, by = 0.015)

# Colon gene expression data
# URL: https://www.kaggle.com/masudur/colon-cancer-gene-expression-data
da <- data.table::fread("data/colon.csv") %>% 
  select(-1) %>% 
  mutate(Class = as.factor(Class))

prop.table(table(da$Class))

# Splitting in train and test data ------------------------------------
set.seed(2019)
split <- da %>% 
  mutate(set = ifelse(
    runif(nrow(.)) > 0.75, "train", "test"))

split %>% 
  count(set) %>% 
  mutate(prop = scales::percent( n/sum(n)))

train <- split %>% filter(set == "train") %>% select(-set)
test <- split %>% filter(set == "test") %>% select(-set)

rf <- randomForest::randomForest(Class ~ ., train)

sum(rf %>% randomForest::importance() > 0)

my_fit_lam <- function(lam){
  m0 <- RRF::RRF(Class ~ ., train,coefReg = lam)
  d <- length(m0$feaSet)  
  pred = predict(m0, test)
  return(list(par = d, tx = sum(pred == test$Class)/nrow(da)))
}



my_fit_gamma <- function(lam){
  gamma_v <- lam
  rm0 <- RRF::RRF(Class ~ ., train)
  impRF <- rm0$importance 
  impRF <- impRF[,"MeanDecreaseGini"] # get the importance score 
  imp <- impRF/(max(impRF)) 
  coefReg <- (1-gamma_v)  + gamma_v*imp   
  grrf <- RRF(Class ~ ., data = train,  
              coefReg = coefReg)
  
  d <- length(grrf$feaSet)  
  return(list(par = d, tx = sum(pred == test$Class)/nrow(da)))
}


lamb_res <- lambdas %>% purrr::map(my_fit_lam)
gamma_res <- lambdas %>% purrr::map(my_fit_gamma)


df <- data.frame(lambdas, 
                 lamb = lamb_res %>% purrr::map_dbl("par"), 
                 gamma = gamma_res %>% purrr::map_dbl("par")) %>% 
  gather("key", "value", -lambdas) %>% 
  mutate(key = forcats::fct_recode(
    key,
    "GRRF: gamma" = "gamma",
    "RRF: lambda" = "lamb"
  ))


p1 <- df %>% 
  ggplot(aes(y = value, lambdas)) +
  facet_wrap(~key, labeller = label_parsed) +
  geom_point() +
  geom_smooth()  + 
  labs(x = 
         expression("Values of the parameters"), 
       y = "Number of variables selected",
       title = "Comparison between the performance of the 
GRRF and RRF when varying their shrinkage parameter") + 
  theme_bw()


p2 <- data.frame(lambdas, 
                 lamb = lamb_res %>% purrr::map_dbl("tx"), 
                 gamma = gamma_res %>% purrr::map_dbl("tx")) %>% 
  gather("key", "value", -lambdas) %>% 
  mutate(key = forcats::fct_recode(
    key,
    "GRRF: gamma" = "gamma",
    "RRF: lambda" = "lamb"
  )) %>% 
  ggplot(aes(y = 1 - value, lambdas)) +
  facet_wrap(~key, labeller = label_parsed) +
  geom_point() +
  geom_smooth(colour = "orange", fill = "white")  + 
  labs(x = 
         expression("Values of the parameters"), 
       y = "Incorrect classification rate",
       caption = 
         "Colon gene expression data,
https://www.kaggle.com/masudur/colon-cancer-gene-expression-data") +
  theme_bw()

library(patchwork)

p1 + p2 + plot_layout(nrow = 2)


# Real gene dataset ---------------------------------------------------
da <-  data.table::fread("data/PV_Data1_ANON.csv") %>% 
  select(-idPhen) %>% select(1:1000)

# Splitting in train and test data ------------------------------------
set.seed(2019)
split <- da %>% 
  mutate(set = ifelse(
    runif(nrow(.)) > 0.75, "train", "test"))

split %>% 
  count(set) %>% 
  mutate(prop = scales::percent( n/sum(n)))

train <- split %>% filter(set == "train") %>% select(-set)
test <- split %>% filter(set == "test") %>% select(-set)


lambdas <- seq(0.01, 0.99, length.out = 50)

msr <- function(model, var){
  res <- var - predict(model, test)
  mean(res^2) 
}


my_fit_lam <- function(lam){
  m0 <- RRF::RRF(log_brd5 ~ ., train, coefReg = lam)
  d <- length(m0$feaSet)  
  return(list(par = d, msr = msr(m0, var = test$log_brd5)))
}



my_fit_gamma <- function(gamma_v){
  rm0 <- RRF::RRF(log_brd5 ~ ., train)
  impRF <- rm0$importance 
  impRF <- impRF[,"IncNodePurity"] # get the importance score 
  imp <- impRF/(max(impRF)) 
  coefReg <- (1 - gamma_v)  + gamma_v*imp   
  grrf <- RRF(log_brd5 ~ ., data = train,  
              coefReg = coefReg)
  
  d <- length(grrf$feaSet)  
  return(list(par = d, msr = msr(grrf, var = test$log_brd5)))
}


lamb_res <- lambdas %>% purrr::map(my_fit_lam)
gamma_res <- lambdas %>% purrr::map(my_fit_gamma)


df <- data.frame(lambdas, 
                 lamb = lamb_res %>% purrr::map_dbl("par"), 
                 gamma = gamma_res %>% purrr::map_dbl("par")) %>% 
  gather("key", "value", -lambdas) %>% 
  mutate(key = forcats::fct_recode(
    key,
    "GRRF: gamma" = "gamma",
    "RRF: lambda" = "lamb"
  ))


p1 <- df %>% 
  ggplot(aes(y = value, lambdas)) +
  facet_wrap(~key, labeller = label_parsed) +
  geom_point() +
  geom_smooth()  + 
  labs(x = 
         expression("Values of the parameters"), 
       y = "Variables selected",
       title = "Comparison between the performance in the test set
of the GRRF and RRF when varying their shrinkage parameter") + 
  theme_bw()


p2 <- data.frame(lambdas, 
           lamb = lamb_res %>% purrr::map_dbl("msr"), 
           gamma = gamma_res %>% purrr::map_dbl("msr")) %>% 
  gather("key", "value", -lambdas) %>% 
  mutate(key = forcats::fct_recode(
    key,
    "GRRF: gamma" = "gamma",
    "RRF: lambda" = "lamb"
  )) %>% 
  ggplot(aes(y = value, lambdas)) +
  facet_wrap(~key, labeller = label_parsed) +
  geom_point() +
  geom_smooth(colour = "orange", fill = "white")  + 
  labs(x = 
         expression("Values of the parameters"), 
       y = "Mean squared residuals",
       caption = 
         "Real gene expression data using the first 1000 variables") +
  theme_bw()


p1 + p2 + plot_layout(nrow = 2)
