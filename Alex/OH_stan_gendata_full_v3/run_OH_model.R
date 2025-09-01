library(rstan)
library(tidybayes)
library(tidyverse)
library(bayesplot)
source("plot_tracker.R")
source("data_gen_v11.R")

# --- Data gen call ---
dataList <- generate_stan_data(
  # sizes & indexing
  N_id = 16,
  N_board = 10,
  N_food = 2,
  seed = 123,
  
  # sampling environment
  total_samples = 120,
  min_board = 8,
  max_board = 14,
  evidence_noise = "multinomial",   # or "deterministic"
  
  # initial OH for test
  mixed_alpha   = 5,                # α for "mixed" group
  uniform_alpha = 0.5,              # α for "uniform" group
  beta_init_mode = "fixed", # or "from_sampling"
  fixed_beta      = c(0.5, 0.5),    # only used if beta_init_mode = "fixed"
  
  # test-phase parameter means (raw scales)
  lambda_raw     = 3,   # log-sensitivity (L = exp(lambda_raw + u_lambda))
  phi_raw        = 0,   # attraction update (P = inv_logit(phi_raw + u_phi))
  stickiness_raw = 1.5, # switch penalty (S = exp(stickiness_raw))
  alpha_last_raw = 1,  # α endpoint (enter squared & negated in sim)
  beta_last_raw  = 1,  # β HV endpoint (enter squared & negated in sim)
  
  # per-ID random effects (λ, α_end, β_end, φ)
  heterogeneity = TRUE,
  sd_lambda   = 0.2,
  sd_alpha_end = 0.30,
  sd_beta_end  = 0.30,
  sd_phi       = 0.50,
  corr_type = "exchangeable",       # "independent", "exchangeable", or "manual"
  rho_exchangeable = 0.2,           # off-diagonal correlation if exchangeable
  corr_matrix_manual = diag(4),     # used only if corr_type = "manual"
  
  # test-phase execution
  max_test_trials = 50,
  break_on_switch = TRUE,
  
  # attractions & eating
  A0 = c(0.8, 0.2),        # initial attractions (HV, LV)
  fixed_A = c(0.8, 0.2),
  p_eat_low  = 0.05,
  p_eat_high = 0.95
)

# Plot test parameters

df <- tibble(T_id = dataList$T_id, treatment = dataList$treatment)
ggplot(df, aes(y=T_id, fill = treatment))+geom_boxplot()+theme_minimal()

plot_tracker(dataList, dataList$trackers$alpha_t, y_lims = c(0,10))
plot_tracker(dataList, dataList$trackers$beta_t, food = 2)
plot_tracker(dataList, dataList$trackers$theta_b11, food = 2)
plot_tracker(dataList, dataList$trackers$theta_b12, food = 2)
plot_tracker(dataList, dataList$trackers$AF_t, food = 1)
plot_tracker(dataList, dataList$trackers$AF_t, food = 2)
plot_tracker(dataList, dataList$trackers$P_S)


####Run Model####
options(mc.cores = parallel::detectCores()-2)

nIter     = 4000
nChains   = 4
nWarmup   = floor(nIter/2)
nThin     = 1

modelFile = 'OH_stan_full_v7.stan'

#Run model
fit = stan(modelFile, 
           data    = dataList, 
           chains  = nChains,
           iter    = nIter,
           warmup  = nWarmup,
           thin    = nThin,
           #init = init_fun,
           #control = list(adapt_delta = 0.95, max_treedepth = 12)
           ) 

####Traceplots####
color_scheme_set("brewer-Spectral")

# 1) Full names (with indices) and base families
full_names <- dimnames(as.array(fit))[[3]]          # e.g. "alpha_s[1]", "beta_s[1,1]", ...
families   <- names(rstan::extract(fit))            # e.g. "alpha_s", "beta_s", ...

# 2) Pick the *first* occurrence for each family
pick_first <- function(par) {
  if (par %in% full_names) return(par)              # scalar
  hits <- full_names[startsWith(full_names, paste0(par, "["))]
  if (length(hits)) hits[1] else NA_character_
}
pars_first <- na.omit(vapply(families, pick_first, character(1)))
pars_first <- unique(c(pars_first, "lp__"))   

# Convert draws to an array (chains x iterations x parameters)
draws_arr <- as.array(fit)
iter <- dim(draws_arr)[2]
n_keep <- 500                                   # <- set what you want
keep <- sort(unique(round(seq(1, iter, length.out = min(n_keep, iter)))))
arr_small <- draws_arr[, keep, , drop = FALSE]

# 2) Trace plots (mixing across chains)
bayesplot::mcmc_trace(arr_small, pars = pars_first)


