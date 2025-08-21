library(rethinking)
library(rstan)
library(tidybayes)
library(tidyverse)
library(reshape2)
source("data_gen_latent_pars8.R")

dataList <- generate_stan_data(N_id = 16 , 
                               mixed_alpha_total = 20, 
                               uniform_alpha_total = 1,
                               fixed_global_counts = NULL,
                               alpha_mode = "decay", #"dirichlet" or "decay"
                               alpha_decay = 0.1, #if alpha_mode = "decay"
                               beta_mode = "global_counts", #"global_counts" or "alpha"
                               rho_a = 1, #if alpha_mode = "dirichlet"
                               rho_b = 1, #if beta_mode = "alpha"
                               lambda_raw = 3,
                               trial_effect = 0,
                               stickiness = 6,
                               break_on_switch = T
)
dataList$T_id
str(dataList)

####Run Model####
options(mc.cores = parallel::detectCores()-2)

nIter     = 2000
nChains   = 4
nWarmup   = floor(nIter/2)
nThin     = 1


modelFile = 'OH_stan_template_v4.stan'

fit = stan(modelFile, 
           data    = dataList, 
           chains  = nChains,
           iter    = nIter,
           warmup  = nWarmup,
           thin    = nThin,
           #init = init_fun
           ) 

fit
#saveRDS(fit,"stanfit_v4.rds")
loadRDS("stanfit_v4.rds")

