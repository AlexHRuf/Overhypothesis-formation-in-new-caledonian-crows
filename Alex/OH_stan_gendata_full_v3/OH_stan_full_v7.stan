//You will need to update your r_stan to versio 2.34 or higher to use the dirichlet_multinomial_lpmf()
//To update to the latest (experimental) version of stan, in your R console run:
//remotes::install_github("stan-dev/rstan@experimental", subdir = "StanHeaders")
//remotes::install_github("stan-dev/rstan@experimental", subdir = "rstan/rstan")

data{
  // Indexing
  int<lower=1> N_test;
  int<lower=1> N_food;             // (=2)
  int<lower=1> N_id;
  int<lower=1> N_board;
  array[N_id] int<lower=1> T_id;   // trials per ID in test
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
  array[N_test] int<lower=1> y_test;              // realised food type
  array[N_test] int<lower=0, upper=1> y_test_ate; // ate? (0/1)
  array[N_test] int<lower=1> id;                  // ID index per test row
  array[N_test] int<lower=1, upper=2> Switch;     // 1=stay(B11), 2=switch(B12)
}

parameters{
  // Preference test - per-ID attractions (exposure)
  real PT_alpha;
  real<lower=0> PT_weight;
  matrix[N_id, N_food] PT_F_value;

  // Sampling overhypotheses (prior for θ per board)
  vector<lower=0>[N_id]      alpha_s;     // scalar α per ID
  array[N_id] simplex[N_food] beta_s;     // mean β per ID

  // Pref + sampling attractions (exposure)
  real ALL_alpha;
  real<lower=0> ALL_weight;
  matrix[N_id, N_food] ALL_F_value;

  // Test phase parameters (raw scale)
  real lambda_raw;             // utility sensitivity
  real phi_raw;                // attraction update rate
  real stickiness_raw;             // Stickiness against switching
  real alpha_last_raw;         // α endpoint (logit-mean for shrink ratio)
  real beta_last_raw;          // β endpoint (logit for HV shrink multiplicative factor)

  // Monotonic schedules (one per ID × trial)
  matrix[N_id, T_max] delta_alpha_raw;  // α schedule logits
  matrix[N_id, T_max] delta_beta_raw;   // β schedule logits

  // Individual random effects for (λ, α-endpoint, β-endpoint, φ)
  matrix[4, N_id] z_ID;                // standard normal
  vector<lower=0>[4] sigma_ID;         // scales
  cholesky_factor_corr[4] Rho_ID;      // correlations
}

transformed parameters{
  // RE offsets: v_ID[j, k], k = 1:λ, 2:α_end, 3:β_end, 4:φ
  matrix[N_id, 4] v_ID = (diag_pre_multiply(sigma_ID, Rho_ID) * z_ID)';

  // Outputs to save
  array[N_test] simplex[2]           Prob_S;       // P(stay,switch) per test row
  array[N_id, T_max] real            a_t;          // α snapshots (scalar)
  array[N_id, N_food, T_max] real    b_t;          // β snapshots (simplex)
  array[N_id, N_food, T_max] real    AF_t;         // attraction snapshots (0..1)

  { // local state (not saved)
    // Per-ID dynamic variables
    array[N_id] real           alpha_t;
    array[N_id] vector[N_food] beta_t;
    array[N_id] vector[N_food] theta_b11;  // pseudo-counts for board 11
    array[N_id] vector[N_food] theta_b12;  // pseudo-counts for board 12

    // Test-phase local counts on B11 (start at 0)
    matrix[N_id, N_food] counts_b11;
    for (j in 1:N_id) counts_b11[j, ] = rep_row_vector(0.0, N_food);

    // Trial counters & padding
    array[N_id] int t_ctr;
    for (j in 1:N_id) {
      t_ctr[j] = 0;
      for (tm in 1:T_max) a_t[j, tm] = 0.0;
      for (f in 1:N_food) for (tm in 1:T_max) b_t[j, f, tm]  = 0.0;
      for (f in 1:N_food) for (tm in 1:T_max) AF_t[j, f, tm] = 0.0;
    }

    // Initial α/β at start of test
    array[N_id] real           a_first;
    array[N_id] vector[N_food] b_first;
    for (j in 1:N_id) {
      a_first[j] = alpha_s[j];
      b_first[j] = beta_s[j];
      alpha_t[j] = a_first[j];
      beta_t[j]  = b_first[j];
      theta_b11[j] = alpha_t[j] * beta_t[j];
      theta_b12[j] = alpha_t[j] * beta_t[j];
    }

    // Attractions initialised from pref+sampling
    matrix[N_id, N_food] ALL_F_run;
    for (j in 1:N_id) for (f in 1:N_food)
      ALL_F_run[j, f] = inv_logit(ALL_F_value[j, f]);

    // Softmax schedules (weights ≥0 sum to 1 over actual trials)
    matrix[N_id, T_max] wA = rep_matrix(0, N_id, T_max);
    matrix[N_id, T_max] wB = rep_matrix(0, N_id, T_max);
    for (j in 1:N_id) {
      vector[T_id[j]] wa = softmax(segment(to_vector(delta_alpha_raw[j]), 1, T_id[j]));
      vector[T_id[j]] wb = softmax(segment(to_vector(delta_beta_raw [j]), 1, T_id[j]));
      wA[j, 1:T_id[j]] = wa';
      wB[j, 1:T_id[j]] = wb';
    }

    // ---------- Test loop ----------
    for (s in 1:N_test) {
      int j = id[s];

      // Current predictive per board (DM predictive via normalisation)
      vector[N_food] P_F_b11 = theta_b11[j] / sum(theta_b11[j]);
      vector[N_food] P_F_b12 = theta_b12[j] / sum(theta_b12[j]);

      // Attractions snapshot
      vector[N_food] A_F;
      for (f in 1:N_food) A_F[f] = ALL_F_run[j, f];

      // Utilities and switch probabilities
      real L = exp(lambda_raw + v_ID[j,1]);
      real S = exp(stickiness_raw);
      real U_b11 = dot_product(P_F_b11, A_F);
      real U_b12 = dot_product(P_F_b12, A_F);

      vector[2] util;
      util[1] = L * U_b11;                 // stay
      util[2] = L * U_b12 - S;    // switch (with penalty)
      Prob_S[s] = softmax(util);

      // Update B11 local counts with realised food
      counts_b11[j, y_test[s]] += 1.0;

      // Advance trial
      t_ctr[j] += 1;

      // Cumulative schedules using PREVIOUS trials (trial 1 stays at start)
      real cumA = (t_ctr[j] >= 2) ? sum(wA[j, 1:(t_ctr[j]-1)]) : 0;
      real cumB = (t_ctr[j] >= 2) ? sum(wB[j, 1:(t_ctr[j]-1)]) : 0;

      // -------- Endpoints for α and β --------
      // α endpoint: shrink-only within (0, a_first]
      //real a_last  = a_first[j] * exp(alpha_last_raw + v_ID[j,2]); //unconstrained -> not working
      //real a_last  = a_first[j] * inv_logit(alpha_last_raw + v_ID[j,2]); // constrained, set prior to normal(0, 1.5) -> centered around 0.5 w. wide tails
      real a_last  = a_first[j] * exp(-square(alpha_last_raw + v_ID[j,2]));  // constrained, set prior to normal(0,1) -> centered near 1 (no change) w, narrorwer tail
      alpha_t[j]   = a_first[j] + (a_last - a_first[j]) * cumA;
      

      // β endpoint: HV can only DECREASE (or stay) multiplicatively from start
      real b_first_HV = b_first[j][1];
      //real p_last_HV  = b_first_HV * exp(beta_last_raw + v_ID[j,3]);  // unconstrained --> not working
      //real p_last_HV  = b_first_HV * inv_logit(beta_last_raw + v_ID[j,3]);  // constrained, set prior to normal(0, 1.5) -> centered around 0.5 w. wide tails
      real p_last_HV  = b_first_HV * exp(-square(beta_last_raw + v_ID[j,3]));  // constrained, set prior to normal(0,1) -> centered near 1 (no change) w. narrower tail
      vector[N_food] b_last;
      b_last[1] = p_last_HV;
      b_last[2] = 1 - p_last_HV;
      beta_t[j] = b_first[j] + (b_last - b_first[j]) * cumB;
      
      // -----------------------------------------------

      // Refresh pseudo-counts (unweighted OH)
      vector[N_food] OH = alpha_t[j] * beta_t[j];
      theta_b11[j] = OH + (counts_b11[j, ]') ; // local counts only on B11
      theta_b12[j] = OH;                       // global OH only on B12

      // Update attractions for chosen food (simple exponential trace)
      real P = inv_logit(phi_raw + v_ID[j,4]);
      for (f in 1:N_food) {
        if (y_test[s] == f) {
          if (y_test_ate[s] + 1 == f) {
            ALL_F_run[j, f] = (1.0 - P) * ALL_F_run[j,f] + P * 0.25;
          } else {
            ALL_F_run[j, f] = (1.0 - P) * ALL_F_run[j,f];
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

model{
  // Preference test (exposure models)
  PT_alpha ~ normal(0, 1);
  PT_weight ~ exponential(1);
  to_vector(PT_F_value) ~ normal(0, 1);

  // Sampling OH priors + collapsed sampling likelihood
  for (j in 1:N_id) {
    alpha_s[j] ~ exponential(1);
    beta_s[j]  ~ dirichlet(rep_vector(1.0, N_food));
  }
  for (j in 1:N_id)
    for (board in 1:N_board)
      target += dirichlet_multinomial_lpmf(
        y_sampling[j, board] | alpha_s[j] * beta_s[j]
      );

  // Random effects (λ, α_end, β_end, φ)
  to_vector(z_ID) ~ normal(0, 1);
  sigma_ID ~ exponential(1);
  Rho_ID ~ lkj_corr_cholesky(4);

  // Pref + sampling attractions (exposure)
  ALL_alpha ~ normal(0, 1);
  ALL_weight ~ exponential(1);
  to_vector(ALL_F_value) ~ normal(0, 1);

  //Schedules and endpoints
  //alpha_last_raw ~ normal(0, 1.5);
  //beta_last_raw  ~ normal(0, 1.5); 
  alpha_last_raw ~ normal(0, 1);
  beta_last_raw  ~ normal(0, 1); 
  to_vector(delta_alpha_raw) ~ normal(0, 1);
  to_vector(delta_beta_raw)  ~ normal(0, 1);

  // Sensitivity / attraction update / stickiness
  lambda_raw ~ normal(0, 1);
  phi_raw    ~ normal(0, 1);
  stickiness_raw ~ normal(0, 1);

  // Pref + sampling totals (exposure models)
  for (j in 1:N_id)
    for (f in 1:N_food) {
      PT_total_eat[j,f]  ~ poisson_log(PT_alpha  + PT_log_dur[j]  + PT_weight  * inv_logit(PT_F_value[j, f]));
      ALL_total_eat[j,f] ~ poisson_log(ALL_alpha + ALL_log_dur[j] + ALL_weight * inv_logit(ALL_F_value[j, f]));
    }

  // Switch choices
  for (s in 1:N_test)
    Switch[s] ~ categorical(Prob_S[s]);
}
