#---------------------------------------------------#
# Regularized Random Forests                        #        
# Bruna Wundervald                                  #        
# April, 2019                                       #        
#---------------------------------------------------#
library(tidyverse) # essential tools
library(RRF) # regularized random forest
library(patchwork)
library(rBayesianOptimization)


rmse <- function(model){
  res <- test$y - predict(model, test)
  sqrt(mean(res^2))
}

da <- read.table("data/fried_added.txt") 

set.seed(2019)

split <- da %>% 
  dplyr::mutate(set = ifelse(
    runif(nrow(.)) > 0.75, "train", "test"))

# Train and test set split 
train <- split %>% filter(set == "train") %>% select(-set)
test <- split %>% filter(set == "test") %>% select(-set)


rf_adj <- function(gamma_v){
  rm0 <- RRF::RRF(y ~ ., train)
  impRF <- rm0$importance 
  impRF <- impRF[,"IncNodePurity"] # get the importance score 
  imp <- impRF/(max(impRF)) #normalize the importance scores into [0,1]
  #gamma_v <- 0.7  
  coefReg <- (1-gamma_v)  + gamma_v*imp   
  grrf <- RRF(y ~ ., data = train, flagReg=1, 
              coefReg = coefReg)
  list(Score = -mean(grrf$mse), Pred = length(grrf$feaSet))
}


lower_bounds <- c(gamma_v = 0)
upper_bounds <- c(gamma_v = 1)
bounds <- list(
  gamma_v = c(lower_bounds[1], upper_bounds[1]))

## Create a grid of values as the input into the BO code
initial_grid <- 
  data.frame(gamma_v = 0.5 + runif(5)/7)

set.seed(2019)
ba_search <- BayesianOptimization(rf_adj,
                                  bounds = bounds,
                                  init_grid_dt = initial_grid, 
                                  init_points = 0, 
                                  n_iter = 10,
                                  acq = "ei", 
                                  kappa = 1, 
                                  eps = 0.0,
                                  verbose = TRUE)


ba_search$Best_Par 
# rm0 <- RRF::RRF(y ~ ., da)
# impRF <- rm0$importance 
# impRF <- impRF[,"MeanDecreaseGini"] # get the importance score 
# imp <- impRF/(max(impRF)) #normalize the importance scores into [0,1]
# gamma_v <- ba_search$Best_Par[1]
# #lambda <- ba_search$Best_Par[1]
# coefReg <- (1-gamma_v) + (gamma_v*imp )
# 
# #ctrl <- trainControl(method = "repeatedcv", repeats = 5)
# 
# 
# grrf <- RRF(Class ~ ., data = da, flagReg = 1,
#             coefReg = coefReg) # 17.74 ): 
# grrf
# 
