#---------------------------------------------------------#
# Regularized Random Forests                              #        
# Simulated data + correlated variables + noisy variables #
# Bruna Wundervald                                        #        
# April, 2019                                             #        
#---------------------------------------------------------#
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
  coefReg <- (1 - gamma)*1 + (gamma*imp)
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


# Applying models ------------------------------------------------------
# models <- da %>% map(run_models)
# saveRDS(models, "code/presentation_version/data/models_corr_true.rds")


#---------------------------------------------------------------
# Calculating the metrics for each model
# -------------------------------------------------------------

models <- readRDS("code/presentation_version/data/models_corr_true.rds")
# Finding number of parameters and RMSE
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

res_pars_formatted <- pars %>% 
  gather("model", "value") %>% 
  mutate(model = forcats::fct_recode(
    model,
    "Using just lambda" = "lambda",
    "Using lambda and gamma" = "gamma",
    "Correlation" = "corre",
    "Using lambda and gamma * correlation" = "comb"
  ))



# Plots! ---------------------------------------------------------------
p1 <- res_formatted %>% 
  filter(!model == "Using just lambda") %>% 
  ggplot(aes(x = value, group = model)) +
  geom_density(aes(fill = model), size = 0.9, 
               colour = "grey35", alpha = 0.7) +
  scale_fill_manual(
    values = c("#919c4c", "#f5c04a", "#e68c7c")) +
  guides(fill = FALSE) +
  xlim(1.5, 3) +
  labs(
    x = "Root mean squared errors", 
    y = "Density",
    fill = 'Variable weighting' 
  )+
  theme_bw() 


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
      expression("Combination of "~gamma~", importance scores, and correlation"),
      "Absolute Correlation",
      expression("GRRF:"~lambda~"= 1 and"~gamma~"= 0.9"))) +
  labs(x = "Number of selected variables",
       y = "Counts",
       fill = 'Variable penalization'
  )  + 
  theme_bw()  +
  theme(legend.position = "bottom",
        legend.direction="vertical")


png("CAIS/img/sim_corr_results.png", 
    width = 200, height = 120, units = 'mm', res = 300)

p1 + p2 + plot_layout(ncol = 1)
dev.off()

# Counting which were the variables selected  --------------------------
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
  filter(!model == "lambda") %>% 
  group_by(model) %>% 
  count(var) %>% 
  arrange(desc(n)) %>% 
  mutate(perc = scales::percent(n/100)) %>% 
  View()

# Table: amount of true variables in each model 
nam <- names(da[[1]]$train)[-6]
ind <- ind %>% 
  mutate(name_var = nam[var])

# Table: correlation of all the variables with the true ones
cors <- cor(da[[1]]$train %>% select(-y))
corrplot::corrplot(cors)

find_cors <- function(var){
  high_cors <- cors[ ,paste(var)][abs(cors[ , paste(var)]) > 0.7]
  true_vars <- paste0("X", 1:5)
  detect <- sum(
    str_detect(names(high_cors), paste0(true_vars, collapse = "|")))
  return(
    ifelse(detect > 0, 1, 0)
  )
}

ind$corr_ind <- ind$name_var %>% map_dbl(find_cors)

write.table(ind, 
            file = "ind_vars_models.txt")

#-----------------------------------------------------------------------