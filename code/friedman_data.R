#---------------------------------------------------#
# Regularized Random Forests                        #        
# Bruna Wundervald                                  #        
# January, 2019                                     #        
#---------------------------------------------------#
library(tidyverse) # essential tools
library(RRF) # regularized random forest
library(rBayesianOptimization) # bayesian optimization 


# Multivariate Adaptive Regression Splines, 1991
# Friedman data simulated in a way that 
# y = beta_1 * sin(pi * x_1 *x_2) + b2 * (x_3 - 0.5)^2 +
# beta_3 * x_4 + beta_4 * x_5 + e; e ~ N(0, 1)

source('https://gist.githubusercontent.com/andrewcparnell/a383e118a27ce809efb09799a3f35939/raw/8201e10249559ab46652aa512a728ba38718aa7e/sim_friedman_simple')

set.seed(2019)
dat = sim_friedman_simple(1000)
fried = data.frame(dat$X, y = dat$y)
head(fried)

# Visualizing ---------------------------------------------------------
fried %>% 
  mutate(X1X2 =  X1 * X2) %>% 
  gather(key = "variable", value = "value", -y) %>% 
  ggplot(aes(y = y, value)) +
  facet_wrap(~variable) + 
  geom_point(colour = "tan3", alpha = 0.8) +
  geom_smooth(color = "grey") +
  theme_bw() +
  labs(x = "Predictor variables", y = "y")
       #title = "Friedman simulated data")

# Introducing noisy variables
set.seed(2019)
A <- matrix(runif(30^2)^2, ncol = 30) 
Sigma <- ((t(A) %*% A)/max(A)) 


sim = MASS::mvrnorm(n = nrow(fried), 
                    mu = seq(0.1, 20, length.out = 30),
                    Sigma = Sigma) %>% 
  as.data.frame()

corrplot::corrplot(cor(sim))

max(cor(sim))
min(cor(sim))

fried <- fried %>% bind_cols(sim)
corrplot::corrplot(cor(fried, method = "spearman"))

# Splitting in train and test data ------------------------------------
set.seed(2019)
split <- fried %>% 
  mutate(set = ifelse(
    runif(nrow(.)) > 0.75, "train", "test"))

split %>% 
  janitor::tabyl(set) %>% 
  mutate(percent = scales::percent(percent))

train <- split %>% filter(set == "train") %>% select(-set)
test <- split %>% filter(set == "test") %>% select(-set)

# Standard Random Forests ----------------------------------------------

(m0 <- randomForest::randomForest(y ~ ., train))

sum(m0$importance > 0)

pred <- predict(m0, newdata = test)
mean((pred - test$y)^2) + 0.8 * sum(m0$importance > 0)


m0 %>% 
  randomForest::importance() %>% 
  as.data.frame() %>% 
  mutate(vars = rownames(.)) %>%
  ggplot(aes(x = reorder(vars, IncNodePurity), IncNodePurity)) +
  geom_linerange(aes(ymin = min(IncNodePurity), ymax = IncNodePurity),
                 position = position_dodge(width = 0.2), size = 1, 
                 colour = 'wheat1') + 
  geom_point(colour = "tan3") + 
  theme_bw() +
  coord_flip() +
  labs(x = "Variables", title = "Importance plot for the Friedman data",
       subtitle = "Standard Random Forest")

# Regularized Random Forests -------------------------------------------
set.seed(2019)
(m1 <- RRF::RRF(y ~ ., train, coefReg = 0.8))

sum(m1$importance > 0)

pred <- predict(m1, newdata = test)
mean((pred - test$y)^2) + 0.8 * sum(m1$importance > 0)

m1 %>% 
  RRF::importance() %>% 
  as.data.frame() %>% 
  mutate(vars = rownames(.)) %>%
  ggplot(aes(x = reorder(vars, IncNodePurity), IncNodePurity)) +
  geom_linerange(aes(ymin = min(IncNodePurity), ymax = IncNodePurity),
                 position = position_dodge(width = 0.2), size = 1, 
                 colour = 'wheat1') + 
  geom_point(colour = "tan3") + 
  theme_bw() +
  coord_flip() +
  labs(x = "Variables", title = "Importance plot for the Friedman data",
       subtitle = "Regularized Random Forest")

# Guided Regularized Random Forests ------------------------------------
set.seed(2019)

impRF <- m1$importance 
impRF <- impRF[,"IncNodePurity"]              # get the importance score 
imp <- impRF/(max(impRF))              # normalize the importance scores

gammas <- seq(0.01, 0.99, length.out = 25)

fit_grrf <- function(gamma){
  coefReg <- (1 - gamma)*1 + (gamma*imp )
  grrf <- RRF(y ~ ., data = train, coefReg = coefReg)
  return(grrf)
}


grrf_grid <- gammas %>% purrr::map(fit_grrf)

importances <- grrf_grid %>% 
  purrr::map(RRF::importance) %>% 
  purrr::map(as.data.frame) %>% 
  purrr::map_df(tibble::rownames_to_column, 'vars', .id = 'name') %>% 
  group_by(vars) %>% 
  summarise(mean_imp = mean(IncNodePurity))


importances %>% 
  as.data.frame() %>% 
  ggplot(aes(x = reorder(vars, mean_imp), mean_imp)) +
  geom_linerange(aes(ymin = min(mean_imp), ymax = mean_imp),
                 position = position_dodge(width = 0.2), size = 1, 
                 colour = 'wheat1') + 
  geom_point(colour = "tan3") + 
  theme_bw() +
  coord_flip() +
  labs(x = "Variables", title = "Mean importance plot for the Friedman data",
       subtitle = "Guided Regularized Random Forest",
       caption = expression(gamma~"varying from 0.01 to 0.99"))


# Selecting the best model 
msr_pars <- function(model){
  vars <- sum(model$importance > 0)
  reg_par <- 0.8
  res <- test$y - predict(model, newdata = test)
  mean(res^2) + (reg_par * vars)
}

msr_mean <- grrf_grid %>% 
  purrr::map_dbl(msr_pars)


# Best lambda
gammas[which.min(msr_mean)]
min(msr_mean)

mod_sel <- grrf_grid[[which.min(msr_mean)]]
sum(mod_sel %>% RRF::importance() > 0)

pred <- predict(mod_sel, newdata = test)
mean((pred - test$y)^2) + 0.8 * sum(mod_sel$importance > 0)

# Guided Regularized Random Forests with Bayesian Optimization ---------

msr_pars <- function(model){
  vars <- sum(model$importance > 0)
  reg_par <- 0.8
  res <- test$y - predict(model, newdata = test)
  mean(res^2) + (reg_par * vars)
}

rf_adj <- function(gamma_v, lambda){
  rm0 <- RRF::RRF(y ~ ., train)
  impRF <- rm0$importance 
  impRF <- impRF[,"IncNodePurity"] 
  imp <- impRF/(max(impRF)) 
  coefReg <- (1-gamma_v)*lambda + gamma_v*imp   
  grrf <- RRF(y ~ ., data = train, coefReg = coefReg)
  list(Score = -msr_pars(grrf), Pred = 0)
}

lower_bounds <- c(gamma_v = 0, lambda = 0)
upper_bounds <- c(gamma_v = 1, lambda = 1)

bounds <- list(
  gamma_v = c(lower_bounds[1], upper_bounds[1]),
  lambda = c(lower_bounds[2], upper_bounds[2]))

## Create a grid of values as the input into the BO code
initial_grid <- 
  data.frame(gamma_v = 0.5 + seq(-0.2, 0.2, by = 0.0335),
             lambda = 0.7 + seq(-0.2, 0.2, by = 0.0335))

set.seed(2019)
# 30 iterations
ba_search <- BayesianOptimization(rf_adj,
                                  bounds = bounds,
                                  init_grid_dt = initial_grid, 
                                  init_points = 0, 
                                  n_iter = 30,
                                  acq = "ucb", 
                                  kappa = 1, 
                                  eps = 0.0,
                                  verbose = TRUE)


# Best parameters found 
ba_search$Best_Par # gamma = 0.9999 and lambda = 0.5740

rm0 <- RRF::RRF(y ~ ., train)
impRF <- rm0$importance 
impRF <- impRF[,"IncNodePurity"] 
imp <- impRF/(max(impRF)) 
gamma_v <- ba_search$Best_Par[1]
lambda <- ba_search$Best_Par[2]
coefReg <- (1-gamma_v)*lambda  + gamma_v*imp  

grrf <- RRF::RRF(y ~ ., data = train, coefReg = coefReg)
grrf
sum(grrf %>% RRF::importance() > 0)
grrf$importance %>% as.data.frame() %>% 
  mutate(var = rownames(grrf$importance)) %>% 
  arrange(desc(IncNodePurity))

pred <- predict(grrf, newdata = test)
mean((pred - test$y)^2) + 0.8 * sum(grrf$importance > 0)
#----------------------------------------------------------------
# Summary of results 

# RF -------------------------------------------------------------
# # Variables = 35
# MSR_test + 0.8 * #vars = 45.34
# MSR_test = 17.34
# RRF -----------------------------------------------------------
# # Variables = 35
# MSR_test + 0.8 * #vars = 45.30
# MSR_test = 17.3
# GRRF: best model with gamma = 0.9 -----------------------------
# # Variables = 17
# MSR_test + 0.8 * #vars = 30.60
# MSR_test = 17
# Bayesian Optimization for gamma and lambda ---------------------
# gamma = 0.8359 and lambda = 0.7137
# # Variables = 18
# MSR_test + 0.8 * #vars = 31.76
# MSR_test = 17.36
#----------------------------------------------------------------

# Saving results
final_models <- list(
  rf = m0, 
  rrf = m1, 
  grrf = mod_sel,
  grrf_bayes = grrf)

saveRDS(final_models, file = "friedman_models.rds")

results_fc <- function(model){
  pred <- predict(model, newdata = test)
  data.frame(
    vars = sum(model$importance > 0),
    msr = mean((pred - test$y)^2)
  ) %>% 
    mutate(var_pars = (vars * 0.8) + msr)
}

results_table <- final_models %>% 
  purrr::map_dfr(results_fc)

write_rds(results_table, "res_fried.rds")
  
