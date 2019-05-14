source('https://gist.githubusercontent.com/andrewcparnell/a383e118a27ce809efb09799a3f35939/raw/8201e10249559ab46652aa512a728ba38718aa7e/sim_friedman_simple')

fc_corr <- function(da){
  var <- sample(names(da)[-6], size = 1)
  var <- da %>% pull(!!sym(var))
  m <- rnorm(1)
  if(m < 0){
    var^(rpois(1, lambda = abs(m)) + rnorm(1, mean = 0.2))
  } 
  else { var * (m + m^2) }
}

simul_data <- function(){
  dat = sim_friedman_simple(1000)
  fried = data.frame(dat$X, y = dat$y)

  vars <- replicate(25, fc_corr(fried)) %>% 
    as.data.frame()
  
  sim <- replicate(n = 30, rnorm(n = nrow(fried), mean = rnorm(1), 
                                 sd = abs(rnorm(1)))) %>% 
    as.data.frame()
  return(
    fried %>% bind_cols(sim) %>% bind_cols(vars) )
  
}

sim_corr <- function(i){
  
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
  
  selected_vars <- paste(nam, collapse = ' + ')
  
  return(list(
    
    train = train,
    test = test,
    selected_vars = selected_vars,
    corr_vars = corr_vars
    
  ))
}

# Simulating 50 datasets ----------------------------------------------
da100corr <- map(1:50, sim_corr)
saveRDS(da100corr, 
        file = "code/presentation_version/data/da100corr_vars.rds")

