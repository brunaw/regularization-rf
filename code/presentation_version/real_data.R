#----------------------------------------------------------#
# Regularized Random Forests                               #     
# Real race horses data                                    #
# Bruna Wundervald                                         #        
# April, 2019                                              #        
#----------------------------------------------------------#
library(tidyverse)
library(RRF)
library(patchwork)

# Real gene data ------------------------------------------------------
da <-  data.table::fread("data/PV_Data1_ANON.csv") %>% 
  select(-idPhen) 

dim(da)
# Number of variables in each model 
mtry = 3582
# Marginal correlations to the response --------------------------------
corr_fc <- function(var){
  cor(da$log_brd5, y = da %>% pull(var))
}

# Splitting in train and test data -------------------------------------
set.seed(2019)
split <- da %>%
  mutate(set = ifelse(
    runif(nrow(.)) < 0.75, "train", "test"))

split %>%
  janitor::tabyl(set) %>%
  mutate(percent = scales::percent(percent))

train <- split %>% filter(set == "train") %>% select(-set)
test <- split %>% filter(set == "test") %>% select(-set)

correlations <- readRDS("code/correlations.rds")

# Applying marginal correlation as the variable weights ----------------
nam <- names(da)[-1]

length(nam[abs(correlations) > 0.15] ) # 3582 variables

selected_vars <- nam[abs(correlations) > 0.15] %>% 
  paste(collapse = ' + ')

corr_weight <- correlations[abs(correlations) > 0.15]

# Split data in different train and test sets --------------------------

sep <- function(i, data, correlations){
  
  split <- da %>% 
    mutate(set = ifelse(
      runif(nrow(.)) < 0.75, "train", "test"))
  
  train <- split %>% filter(set == "train") %>% select(-set)
  test <- split %>% filter(set == "test") %>% select(-set)
  
  selected_vars <- nam[abs(correlations) > 0.10] %>% 
    paste(collapse = ' + ')
  
  corr_weight <- correlations[abs(correlations) > 0.10]
  
  return(list(
    
    train = train,
    test = test,
    selected_vars = selected_vars,
    corr_vars = corr_weight
    
  ))
  
}

# Don't run again
# real_data_50s <- 1:50 %>% map(sep, data = da, correlations = correlations)
# 
# saveRDS(real_data_50, 
#         file = "code/presentation_version/data/real_data_50s.rds")

# real_data_50 <- read_rds(
#  "code/presentation_version/data/real_data_50.rds")

# Splitting into smaller datasets --------------------------------------
real_1_10 <- real_data_50[1:10]
real_11_20 <- real_data_50[11:20]
real_21_30 <- real_data_50[21:30]
real_31_40 <- real_data_50[31:40]
real_41_50 <- real_data_50[41:50]

saveRDS(real_1_10, 
        file = "code/presentation_version/data/real_1_10.rds")
saveRDS(real_11_20, 
        file = "code/presentation_version/data/real_11_20.rds")
saveRDS(real_21_30, 
        file = "code/presentation_version/data/real_21_30.rds")
saveRDS(real_31_40, 
        file = "code/presentation_version/data/real_31_40.rds")
saveRDS(real_41_50, 
        file = "code/presentation_version/data/real_41_50.rds")

# Running models all ---------------------------------------------------

run_models <- function(da, corr_weight = corr_weight, mtry = 3582,
                       selected_vars = selected_vars) {
  train <- da$train
  test <- da$test
  
  rmse <- function(model){
    res <- test$log_brd5 - predict(model, test)
    sqrt(mean(res^2))
  }
  
  form <- paste("log_brd5 ~", selected_vars) %>% as.formula()
  
  # Model 1: fixed lambda 
  ct_lambda <- RRF(form, 
                   data = train, 
                   coefReg =  0.8,
                   mtry = mtry)
  
  # Model 2: guided random forest
  impRF <- ct_lambda$importance
  impRF <- impRF[,"IncNodePurity"]    
  imp <- impRF/(max(impRF))
  
  gamma <-  0.7
  coefReg <- (1 - gamma)*1 + (gamma*imp )
  ct_gamma <- RRF(form, data = train, 
                  coefReg = coefReg, 
                  mtry = mtry)
  
  # Model 3: using the correlation 
  corre <- RRF(form, data = train, 
               coefReg =  corr_weight, 
               mtry = mtry)
  
  # Model 4: combining model
  comb <- data.frame(imp = c(impRF)/(max(c(impRF))),
                     corr = abs(corr_weight)) %>% 
    dplyr::mutate(new_weigth = ifelse(corr > 0.5, 
                                      gamma*imp*corr,
                                      gamma*imp*0.2))
  
  comb_model <- RRF(form, data = train, 
                    coefReg = comb$new_weigth,
                    mtry = mtry)
  
  # Saving only what we need 
  rmse <- list(ct_lambda, ct_gamma, corre, comb_model) %>% 
    map(rmse)
  
  pars <- list(ct_lambda, ct_gamma, corre, comb_model) %>% 
    map(~{ .x$importance})
  
  return(list(
    rmse = rmse,      # RMSE for each model
    pars = pars       # Number of parameters for each model
  ))
}


# Don't run: too slow 
# Next steps are paralellized ------------------------------------------
# Done for all 5 datasets. 
library(furrr)

plan(multiprocess)

results <- future_map(
  1:10, 
  ~{
    models <- run_models(real_1_10[[.x]], selected_vars = selected_vars,
                         corr_weight = corr_weight, mtry = mtry)
    saveRDS(models, 
            paste0(
              "code/presentation_version/data/model_real_dim", .x,  ".rds"))
    remove(models)
    print(.x)
  })
#-----------------------------------------------------------------------
# Obtaining the results ------------------------------------------------
ls_fil <- list.files("code/presentation_version/data/final") 

res <- ls_fil %>%
  paste0("code/presentation_version/data/final/", .) %>% 
  map(read_rds) 

rmse <-   data.frame(
  rmse = res %>% 
    map("rmse") %>% 
    flatten() %>% unlist(),
  model = c("lambda", "gamma", "corre", "comb")
) %>%   
  mutate(model = forcats::fct_recode(
    model,
    "Using just lambda" = "lambda",
    "Using lambda and gamma" = "gamma",
    "Correlation" = "corre",
    "Using lambda and gamma * correlation" = "comb"
  ))

pars <- res %>% 
  map("pars") %>% 
  flatten() %>% 
  map_df(as.data.frame) %>% 
  mutate(model = rep(
    rep(c("lambda", "gamma", "corre", "comb"), each = mtry), 50),
    var = rep(res[[1]]$pars[[1]] %>% rownames(), 4*50),
    model_run = rep(1:50, each = mtry * 4)
  ) %>% 
  mutate(model = forcats::fct_recode(
    model,
    "Using just lambda" = "lambda",
    "Using lambda and gamma" = "gamma",
    "Correlation" = "corre",
    "Using lambda and gamma * correlation" = "comb"
  )) %>% 
  group_by(model_run, model) %>% 
  filter(IncNodePurity > 0) %>% 
  count()



p1 <- rmse %>% 
  filter(!model == "Using just lambda") %>% 
  ggplot(aes(x = rmse, group = model)) +
  geom_density(aes(fill = model), size = 0.9, 
               colour = "grey35", alpha = 0.7) +
  scale_fill_manual(
    values = c("#919c4c", "#f5c04a", "#e68c7c")) +
  guides(fill = FALSE) +
  labs(x = "Root mean squared errors in the test set",
       y = "Density",
       fill = 'Variable weighting') +
  theme_bw() 



p2 <- pars %>% 
  filter(!model == "Using just lambda") %>% 
  ggplot(aes(x = n, group = model)) +
  geom_histogram(aes(fill = model), colour = "grey35",
                 bins = 200,
                 alpha = 0.7, position = "identity") +
  scale_y_continuous(breaks = scales::pretty_breaks())+
  scale_fill_manual(
    guide = guide_legend(),
    values = c("#919c4c", "#f5c04a", "#e68c7c"),
    labels = c(
      expression("Combination of "~gamma~", importance scores, and correlation"),
      "Absolute Correlation",
      expression("GRRF:"~lambda~"= 1 and"~gamma~"= 0.7"))
  ) +
  labs(x = "Number of selected variables", 
       y = "Counts",
       fill = 'Variable penalization')  + 
  xlim(0, 1850) +
  theme_bw()  +
  theme(legend.position = "bottom",
        legend.direction="vertical")


png("CASI/img/real_data_results.png", 
    width = 200, height = 120, units='mm', res = 300)

p1 + p2 + plot_layout(ncol = 1)
dev.off()

# Summary of results
bind_cols(pars %>% 
            group_by(model) %>% 
            summarise(max = max(n), 
                      min = min(n),
                      mean = mean(n)),
          rmse %>% 
            group_by(model) %>% 
            summarise(max_rmse = max(rmse), 
                      min_rmse = min(rmse),
                      mean_rmse = mean(rmse)) %>% 
            select(-model)) %>% 
  write.table("CASI/data/real_res.txt")


# Importance plots -----------------------------------------------------
imp <- res %>% 
  map("pars") %>% 
  flatten() %>% 
  map_df(as.data.frame) %>% 
  mutate(model = rep(
    rep(c("lambda", "gamma", "corre", "comb"), each = mtry), 50),
    var = rep(res[[1]]$pars[[1]] %>% rownames(), 4*50),
    model_run = rep(1:50, each = mtry*4)
  ) %>% 
  group_by(var, model) %>% 
  summarise(imp = mean(IncNodePurity)) %>%
  group_by(model) %>% 
  mutate( var_code = paste0("var_", 1:n())) %>% 
  filter(!model == "lambda") %>% 
  ungroup() %>% 
  mutate_at(vars(model), as.factor)


levels(imp$model) <- c(
  "Combination",
  "Abs.Correlation",
  expression(lambda~"= 1 and"~gamma~"= 0.7"))

png("CASI/img/mean_imp.png", 
    width = 300, height = 120, units='mm', res = 300)

imp %>% 
  group_by(model) %>% 
  arrange(desc(imp)) %>% 
  slice(1:15) %>% 
  ggplot(aes(x = reorder(var_code, imp), y = imp)) +
  geom_linerange(aes(
    ymin = min(imp), ymax = imp, 
    x = reorder(var_code, imp)),
    position = position_dodge(width = 0.2), size = 0.8, 
    colour = '#e68c7c') +
  geom_point(colour = "#e68c7c", size = 2)+
  facet_wrap(~model, scales = "free", labeller = label_parsed) +
  coord_flip() +
  labs(y = "Mean importance measures", 
       x = "Top 15 most important variables") +
  theme_bw(18)

dev.off()  

# Modeling with the best 15 variables of each fit ----------------------
pars_nf <- res %>% 
  map("pars") %>% 
  flatten() %>% 
  map_df(as.data.frame) %>% 
  mutate(model = rep(
    rep(c("lambda", "gamma", "corre", "comb"), each = mtry), 50),
    var = rep(res[[1]]$pars[[1]] %>% rownames(), 4*50),
    model_run = rep(1:50, each = mtry*4)
  ) %>% 
  group_by(model_run, model) %>% 
  filter(IncNodePurity > 0) %>% 
  count() %>% 
  ungroup()

pars_nf %>% 
  select(model_run, model, n) %>% 
  group_by(model_run) %>% 
  spread(key = model, value = n) %>% 
  mutate(win = comb - gamma) %>% 
  ungroup() %>% 
  summarise(s = mean(win))

vars <- imp %>% 
  group_by(model) %>% 
  arrange(desc(imp)) %>% 
  slice(1:15) %>%   
  pull(var)


comb_f <- paste("log_brd5 ~ ", paste0(vars[1:15], collapse = ' + ')) %>% 
  as.formula()
corre_f <- paste("log_brd5 ~ ", paste0(vars[15:30], collapse = ' + '))%>% 
  as.formula()
gamm_f <- paste("log_brd5 ~ ", paste0(vars[30:45], collapse = ' + '))%>% 
  as.formula()



fc <- function(train, test, 
               comb_f = comb_f, 
               gamm_f = gamm_f, 
               corre_f = corre_f){
  rf_g <- randomForest::randomForest(gamm_f, data = train)
  rf_c <- randomForest::randomForest(comb_f, data = train)
  rf_cor <- randomForest::randomForest(corre_f, data = train)
  
  rmse <- function(model, test){
    res <- test$log_brd5 - predict(model, test)
    sqrt(mean(res^2))
  }
  m1 <- rmse(rf_g, test)
  m2 <- rmse(rf_c, test)
  m3 <- rmse(rf_cor, test)
  
  return(list(gamma = m1, comb = m2, corre = m3))
}

# All datasets
da_final <- c(real_1_10, real_11_20, real_21_30, real_31_40, real_41_50)

test <- map2(
  .x = da_final %>% map("train"),
  .y = da_final %>% map("test"), 
  ~{
    fc(train = .x, test = .y, comb_f = comb_f, gamm_f = gamm_f,
       corre_f = corre_f)
  }
)

res_rerun <- data.frame(
  gamma = test %>% map_dbl("gamma"),
  comb = test %>% map_dbl("comb"),
  corre = test %>% map_dbl("corre")
) %>% 
  gather(model, rmse)

# Saving results -------------------------------------------------------
res_rerun %>% 
  group_by(model) %>% 
  summarise(s = mean(rmse), 
            m = max(rmse), 
            min = min(rmse)) %>% 
  write.table("CASI/data/summary_rerun.txt")

# ----------------------------------------------------------------------
