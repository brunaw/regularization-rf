#---------------------------------------------------#
# Regularized Random Forests                        #        
# Bruna Wundervald                                  #        
# April, 2019                                       #        
#---------------------------------------------------#
library(tidyverse) # essential tools
library(RRF) # regularized random forest

# Multivariate Adaptive Regression Splines, 1991
# Friedman data simulated in a way that 
# y = beta_1 * sin(pi * x_1 *x_2) + b2 * (x_3 - 0.5)^2 +
# beta_3 * x_4 + beta_4 * x_5 + e; e ~ N(0, 1)

source('https://gist.githubusercontent.com/andrewcparnell/a383e118a27ce809efb09799a3f35939/raw/8201e10249559ab46652aa512a728ba38718aa7e/sim_friedman_simple')

seed <- function() set.seed(2019)
seed()
dat = sim_friedman_simple(1000)
fried = data.frame(dat$X, y = dat$y)

# Create a few categorical variables
# seed()
# cuts <- runif(n = 5, min = min(dat$y), max = max(dat$y)-5)

# fried <- cuts %>% 
#   purrr::map(
#     ~{
#       data.frame(as.factor(ifelse(fried$y < .x, 0, 1))) %>% 
#         setNames(paste0("var_", round(.x))) 
#       }
#     ) %>% 
#   bind_cols(fried)

# Introducing uncorrelated noisy variables
# seed()

# A <- matrix(runif(50^2)^2, ncol = 50)
# Sigma <- ((t(A) %*% A)/max(A))
# 
# sim = MASS::mvrnorm(n = nrow(fried), 
#                     mu = seq(0.1, 30, length.out = 50),
#                     Sigma = Sigma) %>% 
#   as.data.frame()

seed()
sim <- replicate(n = 95, rnorm(n = nrow(fried), mean = rnorm(1), 
                               sd = abs(rnorm(1)))) %>% 
  as.data.frame()
corrplot::corrplot(cor(sim))

# # Introducing variables that are a function of the predictors
# fc <- function(){
#   nam <- paste0("X", 1:5)
#   ss <- sample(nam, size = 1)
#   power <- rnorm(1)
#   fried %>% pull(!!sym(ss)) %>% .^power
# }  
# 
# var_fc <- replicate(n = 40, fc()) %>% 
#   as.data.frame() %>% 
#   setNames(paste0("V", 51:90))
# 
# fried <- fried %>% bind_cols(sim) %>% bind_cols(var_fc)
fried <- fried %>% bind_cols(sim) 
corrplot::corrplot(cor(fried %>% select_if(is.numeric),
                       method = "spearman"))

da = fried
swrite.table(fried, file = "data/fried_added.txt")

#-----------------------------------------------------------------------