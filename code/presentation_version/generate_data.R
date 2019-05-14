#---------------------------------------------------#
# Regularized Random Forests                        #        
# Bruna Wundervald                                  #        
# April, 2019                                       #        s
#---------------------------------------------------#
library(tidyverse) # essential tools

source('https://gist.githubusercontent.com/andrewcparnell/a383e118a27ce809efb09799a3f35939/raw/8201e10249559ab46652aa512a728ba38718aa7e/sim_friedman_simple')


simul_data <- function(){
  dat = sim_friedman_simple(1000)
  fried = data.frame(dat$X, y = dat$y)
  sim <- replicate(n = 95, rnorm(n = nrow(fried), mean = rnorm(1), 
                                 sd = abs(rnorm(1)))) %>% 
    as.data.frame()
  return(
    fried %>% bind_cols(sim) )
  
}


sim <- function(data){
  
  da <- simul_data()
  nam <- names(da)[-6]
  
  corr_fc <- function(var){
    cor(da$y, y = da %>% pull(var), method = "spearman")
  }
  
  split <- da %>% 
    mutate(set = ifelse(
      runif(nrow(.)) > 0.75, "test", "train"))
  
  # Train and test set split 
  train <- split %>% filter(set == "train") %>% select(-set)
  test <- split %>% filter(set == "test") %>% select(-set)
  
  corr_vars <- nam %>% 
    purrr::map_dbl(corr_fc)
  
  #selected_vars <- nam[abs(corr_vars) > 0.01] %>% 
  #  paste(collapse = ' + ')
  
  selected_vars <- paste(nam, collapse = ' + ')
  
  return(list(
    
    train = train,
    test = test,
    selected_vars = selected_vars,
    corr_vars = corr_vars
    
  ))
}

# Simulating 100 datasets ----------------------------------------------
da100 <- map(1:100, sim)
saveRDS(da100, file = "code/presentation_version/data/da100.rds")

# res %>% 
# ggplot(aes(x = rmse)) +
#   geom_density() +
#   xlim(1.9, 2.7) +
#   theme_bw()
# 
# res %>% 
#   ggplot(aes(x = pars)) +
#   geom_histogram() +
#   theme_bw()
