#---------------------------------------------------#
# Regularized Random Forests                        #        
# Bruna Wundervald                                  #        
# April, 2019                                       #        
#---------------------------------------------------#
library(tidyverse) # essential tools
library(RRF) # regularized random forest
library(patchwork)
#library(rBayesianOptimization)

rmse <- function(model){
  res <- test$y - predict(model, test)
  sqrt(mean(res^2))
}


pars <- function(model){
  length(model$feaSet) 
}

da <- readRDS("code/presentation_version/data/da100.rds")

# da <- da[[1]]

run_models <- function(da) {
  
  train <- da$train
  test <- da$test
  
  rmse <- function(model){
    res <- test$y - predict(model, test)
    sqrt(mean(res^2))
  }
  
  form <- paste("y ~", da$selected_vars) %>% as.formula()
  
  # Model 1: fixed lambda 
  ct_lambda <- RRF(form, 
                   data = train, 
                   coefReg =  0.8,
                   mtry = length(da$corr_vars))
  
  # Model 2: guided random forest
  
  impRF <- ct_lambda$importance
  impRF <- impRF[,"IncNodePurity"]    
  imp <- impRF/(max(impRF))
  
  gamma = 0.9
  coefReg <- (1 - gamma)*1 + (gamma*imp )
  ct_gamma <- RRF(form, data = train, coefReg = coefReg, 
                  mtry = length(da$corr_vars))
  
  # Model 3: using the correlation 
  corre <- RRF(form, data = train, 
               coefReg =  abs(da$corr_vars), 
               mtry = length(da$corr_vars))
  # Model 4: combining model
  
  comb <- data.frame(imp = c(impRF)/(max(c(impRF))),
                  corr = abs(da$corr_vars)) %>% 
    dplyr::mutate(new_weigth = ifelse(corr > 0.5, 
                                      gamma*imp*corr,
                                      gamma*imp*0.3))
  
  
 
  comb_model <- RRF(form, data = train, 
                    coefReg = comb$new_weigth,
                    mtry = length(da$corr_vars))
  
  return(list(
    ct_lambda = ct_lambda, 
    ct_gamma = ct_gamma, 
    corre = corre, 
    comb_model = comb_model
  ))
}


# Applying models

models <- da %>% map(run_models)
saveRDS(models, "code/presentation_version/data/models.rds")

#---------------------------------------------------------------
# Calculating the metrics for each model
# -------------------------------------------------------------

# Finding n-parameters 
pars <- models %>% 
  map_df(~{
    data.frame(
      lambda = length(.x$ct_lambda$feaSet),
      gamma = length(.x$ct_gamma$feaSet),
      corre = length(.x$corre$feaSet),
      comb = length(.x$comb_model$feaSet))
  })


# Finding mse
mse <- map2_df(
  .x = models,
  .y = da %>% map("test"),
  ~{
    data.frame(
      lambda = sqrt(mean((.y$y - predict(.x$ct_lambda, .y))^2 )),
      gamma = sqrt(mean((.y$y - predict(.x$ct_gamma, .y))^2 )),
      corre = sqrt(mean((.y$y - predict(.x$corre, .y))^2 )),
      comb = sqrt(mean((.y$y - predict(.x$comb_model, .y))^2 ))
    )
  }
)

res_formatted <- mse %>% 
  gather("model", "value") %>% 
  mutate(model = forcats::fct_recode(
    model,
    "Using just lambda" = "lambda",
    "Using lambda and gamma" = "gamma",
    "Correlation" = "corre",
    "Using lambda and gamma * correlation" = "comb"
  ))


p1 <- res_formatted %>% 
  filter(!model == "Using just lambda") %>% 
  ggplot(aes(x = value, group = model)) +
  geom_density(aes(fill = model), size = 0.9, 
               colour = "grey35", alpha = 0.7) +
  scale_fill_manual(
    values = c("#919c4c", "#f5c04a", "#e68c7c")) +
  guides(fill = FALSE) +
  labs(x = "Root mean squared errors", y = "Density",
       fill = 'Variable weighting', 
       title = "Comparing different weightings for the GRRF",
       subtitle = "Results obtained with 100 different simulated datasets") +
  theme_bw() 


res_pars_formatted <- pars %>% 
  gather("model", "value") %>% 
  mutate(model = forcats::fct_recode(
    model,
    "Using just lambda" = "lambda",
    "Using lambda and gamma" = "gamma",
    "Correlation" = "corre",
    "Using lambda and gamma * correlation" = "comb"
  ))


p2 <- res_pars_formatted %>% 
  filter(!model == "Using just lambda") %>% 
  ggplot(aes(x = value, group = model)) +
  geom_histogram(aes(fill = model), colour = "grey35",
                 bins = 35,
                 alpha = 0.7, position = "identity") +
  scale_fill_manual(
    guide = guide_legend(),
    values = c("#919c4c", "#f5c04a", "#e68c7c"),
    labels = c(
      # expression(lambda~"= 0.8"),
      expression("GRRF: combination of "~lambda~", "~gamma~" and correlation"), 
      "Correlation as reg. parameter",
      expression("GRRF:"~lambda~"= 0.8 and"~gamma~"= 0.9"))
      #  expression("GRRF: combination of "~lambda~", "~gamma~" and #correlation"))) +
  )+
  labs(x = "Number of selected variables", y = "Counts",
       fill = 'Variable weighting', 
       caption = 
         "Friedman data")  + 
  theme_bw()  +
  theme(legend.position = "bottom",
        legend.direction="vertical")


p1 + p2 + plot_layout(ncol = 1)

# Count which were the variables selected 
ind <- 
  data.frame(
    var = c(
      models %>% 
        map(~{ .x$ct_gamma$feaSet } ) %>% unlist(),
      models %>% 
        map(~{ .x$corre$feaSet } ) %>% unlist(),
      models %>% 
        map(~{ .x$comb_model$feaSet } ) %>% unlist()),
    model = c(
      rep("gamma", sum(pars$gamma)),
      rep("corre", sum(pars$corre)),
      rep("comb", sum(pars$comb)))
  )


ind %>% 
  group_by(model) %>% 
  count(var) %>% 
  arrange(desc(n)) %>% 
  mutate(perc = scales::percent(n/100)) %>% 
  view()

write.table(ind, 
            file = "code/presentation_version/data/ind_vars_models.txt")
  


