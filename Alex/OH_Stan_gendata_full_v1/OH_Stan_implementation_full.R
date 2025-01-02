library(rethinking)
library(rstan)
library(tidybayes)
library(tidyverse)
library(reshape2)
source("gen_srtan_data.R")

dataList <- generate_stan_data(N_id = 6)
str(dataList)

####Run Model####
options(mc.cores = parallel::detectCores()-2)

nIter     = 2000
nChains   = 1
nWarmup   = floor(nIter/2)
nThin     = 1


modelFile = 'OH_stan_template_full_alexis.stan'

fit = stan(modelFile, 
           data    = dataList, 
           chains  = nChains,
           iter    = nIter,
           warmup  = nWarmup,
           thin    = nThin,
           #init = init_fun
           ) 

fit




