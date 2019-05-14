#--------------------------------------------------------#
# Regularized Random Forests                             #    
# Random forests in the simulated data                   #
# Bruna Wundervald                                       #        
# April, 2019                                            #        
#--------------------------------------------------------#
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

da <- readRDS("code/presentation_version/data/da100corr_vars.rds")

# Models for all the variables -----------------------------------------
run_rf <- function(da) {
  
  train <- da$train
  test <- da$test
  rmse <- function(model){
    res <- test$y - predict(model, test)
    sqrt(mean(res^2))
  }
  
  # Standard Random Forest model
  rf <- randomForest::randomForest(y ~ ., 
                                   data = train)

  impRF <- data.frame(rf$importance) %>% 
    mutate(vars = rownames(.))
  return(impRF)
}


models_rf <- da %>% map(run_rf)
saveRDS(models_rf, "code/presentation_version/data/rf_true.rds")

# Models without the correlated variables ------------------------------
run_rf_dim <- function(da) {
  
  train <- da$train
  test <- da$test
  vars <- names(train)[1:37][-6]
  form <- paste0("y ~", 
                 paste(vars, collapse = " + ")) %>% 
    as.formula()
  rmse <- function(model){
    res <- test$y - predict(model, test)
    sqrt(mean(res^2))
  }
  
  
  # RF
  rf <- randomForest::randomForest(form, 
                                   data = train)
  
  impRF <- data.frame(rf$importance) %>% 
    mutate(vars = rownames(.))
  
  return(impRF)
}


# models_rf_dim <- da %>% map(run_rf_dim)
# saveRDS(models_rf_dim, "code/presentation_version/data/rf_true_dim.rds")
models_rf <- readRDS("code/presentation_version/data/rf_true.rds")
models_rf_dim <- readRDS("code/presentation_version/data/rf_true_dim.rds")

# Evaluating results ---------------------------------------------------
lens <- data.frame(
  full = models_rf %>% 
    map_dbl(nrow),
  dim = models_rf_dim %>% 
    map_dbl(nrow))

var_imps <- bind_rows(
  models_rf %>% 
  map_df(data.frame) %>% 
  mutate(model_prov = rep(1:50, each = 60),
         rf = "Including the correlated variables"),
  
  models_rf_dim %>% 
    map_df(data.frame) %>% 
    mutate(model_prov = rep(1:50, each = 36),
           rf = "Without the correlated variables")
  
  ) %>% 
  group_by(rf, vars) %>% 
  summarise(imp = mean(IncNodePurity)) %>% 
  group_by(rf) %>% 
  arrange(desc(imp)) %>% 
  slice(1:25) %>% 
  ungroup() %>% 
  mutate_at(vars(rf), as.factor) %>% 
  mutate(rf = forcats::fct_relevel(
    rf,
    c("Without the correlated variables",
      "Including the correlated variables")
  ))

glimpse(var_imps)  

# Plotting the comparison ----------------------------------------------
png("CASI/img/rf_comparison.png", 
    width = 300, height = 300, units='mm', res = 300)

var_imps %>% 
  ggplot(aes(x = reorder(vars, imp), y = imp)) +
  geom_linerange(aes(
    ymin = min(imp), ymax = imp, 
    x = reorder(vars, imp)),
    position = position_dodge(width = 0.2), size = 0.8, 
    colour = '#e68c7c') +
  geom_point(colour = "#e68c7c", size = 2)+
  facet_wrap(~rf, scales = "free") +
  coord_flip() +
  labs(y = "Importance measures", 
       x = "Top 25 most important variables") +
  theme_bw(16)
dev.off()

# ---------------------------------------------------------------------