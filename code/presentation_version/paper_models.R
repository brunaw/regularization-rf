#---------------------------------------------------#
# Regularized Random Forests                        #        
# Bruna Wundervald                                  #        
# Arpil, 2019                                       #        
#---------------------------------------------------#
library(tidyverse) # essential tools
library(RRF) # regularized random forest
library(patchwork)

da <- read.table("data/fried_added.txt") 

set.seed(2019)

split <- da %>% 
  mutate(set = ifelse(
    runif(nrow(.)) > 0.75, "train", "test"))

# Train and test set split 
train <- split %>% filter(set == "train") %>% select(-set)
test <- split %>% filter(set == "test") %>% select(-set)


# MSE function
rmse <- function(model, var){
  res <- var - predict(model, test)
  sqrt(mean(res^2))
}


# Defining many lambdas for the regularisation
lambdas <- c(0.0001, seq(0.00001, 0.99, length.out = 50))
# RRF
my_fit_lam <- function(lam){
  m0 <- RRF::RRF(y ~ ., data = train, coefReg = lam)
  d <- length(m0$feaSet)  
  return(list(par = d, rmse = rmse(m0, var = test$y)))
}

# GRRF
my_fit_gamma <- function(gamma_v){
  rm0 <- RRF::RRF(y ~ ., train)
  impRF <- rm0$importance 
  impRF <- impRF[,"IncNodePurity"] # get the importance score 
  imp <- impRF/(max(impRF)) 
  coefReg <- (1 - gamma_v)  + gamma_v*imp   
  grrf <- RRF(y ~ ., data = train,  
              coefReg = coefReg)
  
  d <- length(grrf$feaSet)  
  return(list(par = d, rmse = rmse(grrf, var = test$y)))
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
  geom_hline(yintercept = 5, colour = "red") +
  labs(x = 
         expression("Values of the parameters"), 
       y = "Variables selected",
       title = "Comparison between the performance in the test set
of the GRRF and RRF when varying their shrinkage parameter") + 
  theme_bw()


p2 <- data.frame(lambdas, 
                 lamb = lamb_res %>% purrr::map_dbl("rmse"), 
                 gamma = gamma_res %>% purrr::map_dbl("rmse")) %>% 
  gather("key", "value", -lambdas) %>% 
  mutate(key = forcats::fct_recode(
    key,
    "GRRF: gamma" = "gamma",
    "RRF: lambda" = "lamb"
  )) %>% 
  ggplot(aes(y = value, lambdas)) +
  facet_wrap(~key, labeller = label_parsed) +
  geom_point() +
  geom_smooth(colour = "orange", alpha = 0)  + 
  labs(x = 
         expression("Values of the parameters"), 
       y = "RMSE",
       caption = 
         "Friedman data with added noisy variables") +
  theme_bw()


p1 + p2 + plot_layout(nrow = 2)
