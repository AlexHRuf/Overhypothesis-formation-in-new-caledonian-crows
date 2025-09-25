// Combo model: M1 RE & decays + M2 attraction updating + FIXED stickiness
// Requires rstan >= 2.34 for dirichlet_multinomial_lpmf

data {
  // Indexing
  int<lower=1> N_test;
  int<lower=1> N_food;               // (=2)
  int<lower=1> N_id;
  int<lower=1> N_board;
  array[N_id] int<lower=1> T_id;     // trials per ID in test (not strictly needed for decays)
  int<lower=1> T_max;

  // Preference test (exposure model)
  array[N_id, N_food] int PT_total_eat;
  array[N_id] real PT_log_dur;

  // Sampling phase (boards 1..N_board)
  array[N_id, N_board, N_food] int y_sampling;

  // Pref + sampling totals (exposure model)
  array[N_id, N_food] int ALL_total_eat;
  array[N_id] real ALL_log_dur;

  // Test phase sequence
  array[N_test] int<lower=1> y_test;              // realised food type (1..N_food)
  array[N_test] int<lower=0, upper=1> y_test_ate; // ate? (0/1)
  array[N_test] int<lower=1> id;                  // ID index per test row
  array[N_test] int<lower=1, upper=2> Switch;     // 1=stay(B11), 2=switch(B12)
}

parameters {
  // Preference test - per-ID attractions (exposure)
  real PT_alpha;
  real<lower=0> PT_weight;
  matrix[N_id, N_food] PT_F_value;

  // Sampling overhypotheses (prior for θ per board), collapsed in likelihood
  vector<lower=0>[N_id] alpha_s;                  // scalar α per ID (natural scale)
  array[N_id] simplex[N_food] beta_s;             // mean β per ID (natural scale)

  // Pref + sampling attractions (exposure)
  real ALL_alpha;
  real<lower=0> ALL_weight;
  matrix[N_id, N_food] ALL_F_value;

  // Test phase (population means, raw scales)
  real lambda_raw;                                 // utility sensitivity (mean, log-scale)
  real phi_raw;                                    // attraction update rate (mean, logit)

  // Simple exponential decays (Model 1 logic)
  real<lower=0, upper=0.2> d_alpha;               // per-trial decay on α_t
  real<lower=0, upper=0.2> d_beta;                // per-trial decay on β_t[1] (HV share)

  // Individual random effects for (λ, φ) — 2-variate
  matrix[2, N_id] z_ID;                           // standard normal
  vector<lower=0>[2] sigma_ID;                    // scales
  cholesky_factor_corr[2] Rho_ID;                 // correlations
}

transformed parameters {
  // RE offsets: v_ID[j,1]=λ offset, v_ID[j,2]=φ offset
  matrix[N_id, 2] v_ID = (diag_pre_multiply(sigma_ID, Rho_ID) * z_ID)';

  // Saved outputs
  array[N_test] simplex[2]         Prob_S;        // P(stay,switch)
  array[N_id, T_max] real          a_t;           // α snapshots (scalar)
  array[N_id, N_food, T_max] real  b_t;           // β snapshots (simplex)
  array[N_id, N_food, T_max] real  AF_t;          // attractions (0..1)

  { // local (not saved)
    // Dynamic state per ID
    array[N_id] real           alpha_t;           // scalar
    array[N_id] vector[N_food] beta_t;            // simplex
    array[N_id] vector[N_food] theta_b11;         // pseudo-counts board 11
    array[N_id] vector[N_food] theta_b12;         // pseudo-counts board 12

    // Test-phase local counts on B11 (start at 0)
    matrix[N_id, N_food] counts_b11;
    for (j in 1:N_id) counts_b11[j, ] = rep_row_vector(0.0, N_food);

    // Trial counters & padding
    array[N_id] int t_ctr;
    for (j in 1:N_id) {
      t_ctr[j] = 0;
      for (tm in 1:T_max) a_t[j, tm] = 0.0;
      for (f in 1:N_food) for (tm in 1:T_max) {
        b_t[j, f, tm]  = 0.0;
        AF_t[j, f, tm] = 0.0;
      }
    }

    // Initial α/β at start of test from OH prior
    for (j in 1:N_id) {
      alpha_t[j]  = alpha_s[j];
      beta_t[j]   = beta_s[j];
      theta_b11[j] = alpha_t[j] * beta_t[j];
      theta_b12[j] = alpha_t[j] * beta_t[j];
    }

    // Attractions initialised from pref+sampling
    matrix[N_id, N_food] ALL_F_run;
    for (j in 1:N_id) for (f in 1:N_food)
      ALL_F_run[j, f] = inv_logit(ALL_F_value[j, f]);

    // ------------ Test loop ------------
    for (s in 1:N_test) {
      int j = id[s];

      // Predictive per board (normalize pseudo-counts)
      vector[N_food] P_F_b11 = theta_b11[j] / sum(theta_b11[j]);
      vector[N_food] P_F_b12 = theta_b12[j] / sum(theta_b12[j]);

      // Attractions snapshot
      vector[N_food] A_F;
      for (f in 1:N_food) A_F[f] = ALL_F_run[j, f];

      // Utilities and stay/switch probabilities
      real L = exp(lambda_raw + v_ID[j, 1]);  // per-ID sensitivity
      real S = exp(1.5);                      // FIXED stickiness: kappa_raw = 1.5
      real U_b11 = dot_product(P_F_b11, A_F);
      real U_b12 = dot_product(P_F_b12, A_F);

      vector[2] util;
      util[1] = L * U_b11;          // stay B11
      util[2] = L * U_b12 - S;      // switch B12 with fixed penalty
      Prob_S[s] = softmax(util);

      // Update B11 local counts with realised food
      counts_b11[j, y_test[s]] += 1.0;

      // Advance trial and apply simple decays (Model 1 logic)
      t_ctr[j] += 1;
      alpha_t[j]     = (1.0 - d_alpha) * alpha_t[j];     // α decay
      beta_t[j][1]   = (1.0 - d_beta)  * beta_t[j][1];   // β decay on HV share
      beta_t[j][2]   = 1.0 - beta_t[j][1];

      // Refresh pseudo-counts (unweighted OH)
      vector[N_food] OH = alpha_t[j] * beta_t[j];
      theta_b11[j] = OH + (counts_b11[j, ]') ;   // local counts on B11
      theta_b12[j] = OH;                         // OH only on B12

      // Attraction updating (Model 2 style), contingent on choice & eating
      real P = inv_logit(phi_raw + v_ID[j, 2]);  // per-ID RE on φ
      for (f in 1:N_food) {
        if (y_test[s] == f) {
          if (y_test_ate[s] + 1 == f) {
            ALL_F_run[j, f] = (1.0 - P) * ALL_F_run[j, f] + P * 0.25;
          } else {
            ALL_F_run[j, f] = (1.0 - P) * ALL_F_run[j, f];
          }
        }
      }

      // Save snapshots
      a_t[j, t_ctr[j]] = alpha_t[j];
      for (f in 1:N_food) {
        b_t[j, f, t_ctr[j]]  = beta_t[j][f];
        AF_t[j, f, t_ctr[j]] = ALL_F_run[j, f];
      }
    }
  }
}

model {
  // Preference test (exposure models)
  PT_alpha ~ normal(0, 1);
  PT_weight ~ exponential(1);
  to_vector(PT_F_value) ~ normal(0, 1);

  // Sampling OH priors + collapsed sampling likelihood (DM)
  alpha_s ~ exponential(1);
  for (j in 1:N_id) beta_s[j] ~ dirichlet(rep_vector(1.0, N_food));
  for (j in 1:N_id)
    for (board in 1:N_board)
      target += dirichlet_multinomial_lpmf(
        y_sampling[j, board] | alpha_s[j] * beta_s[j]
      );

  // 2D random effects for (λ, φ)
  to_vector(z_ID) ~ normal(0, 1);
  sigma_ID ~ exponential(1);
  Rho_ID ~ lkj_corr_cholesky(2);

  // Pref + sampling attractions (exposure)
  ALL_alpha ~ normal(0, 1);
  ALL_weight ~ exponential(1);
  to_vector(ALL_F_value) ~ normal(0, 1);

  // Decay parameters (Model 1 priors)
  d_alpha ~ beta(2, 330);  // ~0.006 mean; ~75% remains after ~50 trials
  d_beta  ~ beta(2, 330);

  // Population means (raw) for λ and φ
  lambda_raw ~ normal(0, 1);
  phi_raw    ~ normal(0, 1);

  // Pref + sampling totals (exposure models)
  for (j in 1:N_id)
    for (f in 1:N_food) {
      PT_total_eat[j, f]  ~ poisson_log(PT_alpha  + PT_log_dur[j]  + PT_weight  * inv_logit(PT_F_value[j, f]));
      ALL_total_eat[j, f] ~ poisson_log(ALL_alpha + ALL_log_dur[j] + ALL_weight * inv_logit(ALL_F_value[j, f]));
    }

  // Switch choices (uses Prob_S built with fixed S)
  for (s in 1:N_test)
    Switch[s] ~ categorical(Prob_S[s]);
}
