library(rstan); library(tidybayes); library(tidyverse); library(bayesplot)

source("plot_tracker.R")
source("data_gen_simple.R") #stripped down data gen without parameter inputs.
source("data_gen_full.R") #parameter inputs

dataList <- generate_stan_data_full(
  N_id = 20, N_board = 10, seed = 1234,
  total_samples = 120, min_board = 8, max_board = 14,
  
  alpha_mixed = 5,  
  alpha_uniform = 0.5,
  beta_mixed  = c(0.5, 0.5),
  beta_uniform= c(0.5, 0.5),
  
  # RAW means
  lambda_raw_mixed   = 2.8,   # log L
  lambda_raw_uniform = 3.2,
  phi_raw_mixed      = 0.0,   # logit P
  phi_raw_uniform    = 0.0,
  kappa_raw_mixed    = 1.5,   # log K
  kappa_raw_uniform  = 1.5,
  
  heterogeneity = TRUE,
  sd_lambda = 0.1, sd_phi = 0.1, sd_kappa = 0.1,
  corr_type = "exchangeable", rho_exchangeable = 0.2,
  
  max_test_trials = 50, break_on_switch = TRUE,
  A0 = c(0.8, 0.2),
  fixed_A = c(0.8, 0.2),                    # let φ matter
  d_alpha = 0.01, d_beta = 0.01
)
plot_tracker(dataList, dataList$trackers$P_S)

dataList$T_id
df <- tibble(T_id = dataList$T_id, treatment = dataList$treatment)
ggplot(df, aes(y=T_id, fill = treatment))+geom_boxplot()

plot_tracker(dataList, dataList$trackers$alpha_t, y_lims = c(0,10))
plot_tracker(dataList, dataList$trackers$beta_t, food = 2)
plot_tracker(dataList, dataList$trackers$theta_b11, food = 2)
plot_tracker(dataList, dataList$trackers$theta_b12, food = 2)
plot_tracker(dataList, dataList$trackers$AF_t, food = 1)
plot_tracker(dataList, dataList$trackers$AF_t, food = 2)
plot_tracker(dataList, dataList$trackers$P_S)


####Run Model####
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()-2)

nIter     = 4000
nChains   = 4
nWarmup   = floor(nIter/2)
nThin     = 1

modelFile_no_kappa = 'OH_stan_full_v8_no_kappa.stan'
modelFile_global_kappa = 'OH_stan_full_v8_global_kappa.stan'
modelFile_id_kappa = 'OH_stan_full_v8_id_kappa.stan'
modelFile_fixed_kappa = 'OH_stan_full_v8_fixed_kappa.stan'

#Run models
fit_no_kappa = stan(modelFile_no_kappa, data = dataList, chains = nChains, iter = nIter, warmup = nWarmup, thin = nThin) 
fit_global_kappa = stan(modelFile_nostick, data = dataList, chains = nChains, iter = nIter, warmup = nWarmup, thin = nThin) 
fit_id_kappa = stan(modelFile_id_kappa, data = dataList, chains = nChains, iter = nIter, warmup = nWarmup, thin = nThin) 
fit_fixed_kappa = stan(modelFile_fixed_kappa, data = dataList, chains = nChains, iter = nIter, warmup = nWarmup, thin = nThin) 

####Traceplots####
fit <- fit_fixed_kappa

source("traceplot_first.r")
traceplot_first(fit, n_keep = 500, color_scheme = "brewer-Spectral")
