data{
  // General indexing 
  int<lower=1> N_test;        // Total test-phase observations 
  int<lower=1> N_food;        // Total food items (2)
  int<lower=1> N_id;          // Total individuals (16?)
  int<lower=1> N_board;       // Total boards in sampling phase (10)
  int<lower=0> T_id[N_id];    // Trials per crow in test (can be 0 for some)
  int<lower=1> T_max;         // max(T_id) used to pad arrays to rectangular shape

  // Preference test - we can use these to show a priori to manipulation (sampling phase) crows differed in how much the liked the foods on offer
  array[N_id, N_food] int PT_total_eat; // Total per crow high- and low-quality food counts 
  array[N_id] real PT_log_dur;          // Total per crow log duration of pref test; exposure parameter accounting for differences in eating opportunity

  // Sampling phase 
  array[N_id, N_board, N_food] int y_sampling; // Total per crow, per board high- and low-quality food counts (1 = high, 2 = low)

  // Preference test + sampling phase totals 
  array[N_id, N_food] int ALL_total_eat; // Total per crow high- and low-quality food counts across all eating
  array[N_id] real ALL_log_dur;          // Total per crow log duration across all eating; exposure parameter accounting for differences in eating opportunity

  // Test phase 
  array[N_test] int<lower=1> y_test;              // Realised food from board 11 - until switch, will be entirely low-quality (i.e., 2s)
  array[N_test] int<lower=0, upper=1> y_test_ate; // Ate that item? (0/1) Need these for food-item attraction updating
  array[N_test] int<lower=1> id;                  // Individual index
  array[N_test] int<lower=1> trial;               // Trial number
  array[N_test] int<lower=1, upper=2> Switch;     // Switch? 1 = stay on Board 11, 2 = switch to Board 12
}

parameters{
  // Preference test - we estimate this individually, to show reviewers that crows will have a 'fav' food-item a priori to sampling phase manipulation
  real PT_alpha;                   // Intercept i.e., rate of random choice 
  real<lower=0> PT_weight;         // Weight for food-item attraction values i.e., allow attractions to differentially direct eating, because can't assume 'like' or 'dislike' has same effect on food-choice across crows
  matrix[N_id, N_food] PT_F_value; // High and Low quality individual-level attraction value estimates for all food eaten 

  // Sampling overhypotheses 
  vector<lower=0>[N_id] S_gamma;          // Over-overhypothesis
  vector<lower=0>[N_food] S_alpha[N_id];  // Overhypothesis about food item variability
  simplex[N_food] S_beta[N_id];           // Overhypothesis about food item distribution
  simplex[N_food] S_theta[N_id, N_board]; // Local (per board) knowledge about food item variability & distribution

  // Pref + sampling attractions 
  real ALL_alpha;                   // Intercept i.e., rate of random choice 
  real<lower=0> ALL_weight;         // Weight for food-item attraction values i.e., allow attractions to differentially direct eating, because can't assume 'like' or 'dislike' has same effect on food-choice across crows
  matrix[N_id, N_food] ALL_F_value; // High and Low quality individual-level attraction value estimates for all food eaten 

  // Test phase target parameters on raw scale
  real rho_alpha; // OH alpha updating rate (scaler) 
  real rho_beta;  // OH beta updating rate (scaler)
  real lambda;    // Sensitivity rate (weight) to differences in board utility (i.e., prob food informed by OH * attraction to food)
  real phi;       // Food attraction values updating rate (scaler)
  
  // Test phase covariate
  real trial_effect; // Temporal bias toward switching

  // Individual random effects 
  matrix[4, N_id] z_ID;           // Standard normal per target parameter per crow
  vector<lower=0>[4] sigma_ID;    // Standard deviation per target parameter
  cholesky_factor_corr[4] Rho_ID; // Correlations per target parameters
}

transformed parameters{
  // Random effect offsets
  matrix[N_id, 4] v_ID = (diag_pre_multiply(sigma_ID, Rho_ID) * z_ID)';

  // Outputs saved in global
  simplex[2] Prob_S[N_test];           // Stay vs switch by trial -> goes to likelihood
  array[N_id, N_food, T_max] real a_s; // Padded snapshot OH alpha_s (i.e., zeros in unused slots); facilitates array indexing b/c otherwise would be ragged structure
  array[N_id, N_food, T_max] real a_b; // Padded snapshot OH beta_s  (i.e., zeros in unused slots); facilitates array indexing b/c otherwise would be ragged structure

  { // Local scope - nothing below here saved
  
    // Dynamic variable-updates per crow in-test
    array[N_id] vector[N_food] alpha_s;   // OH variability (pseudo-prior counts):
                                          // Dirichlet-style "strengths" representing
                                          // prior experience of each food type.

    array[N_id] vector[N_food] beta_s;    // OH distribution (pseudo-prior probabilities):
                                          // Normalized version of alpha_s (i.e. beta = alpha / sum(alpha)),
                                          // yielding simplex values of belief (0-1, summing to 1).
                                          // Matches Dirichlet logic where alpha defines pseudo-counts
                                          // and beta is the implied belief distribution.

    array[N_id] vector[N_food] theta_b11; // Board 11 belief (pseudo-prior counts):
                                          // Combination of OH (alpha_s * beta_s) with direct counts
                                          // from observed samples on Board 11.

    array[N_id] vector[N_food] theta_b12; // Board 12 belief (pseudo-prior counts):
                                          // Based only on OH (alpha_s * beta_s), without direct
                                          // Board 11 evidence. Represents the 'pure' prior belief
                                          // carried into test for the alternative board.
    
    // Initialise food counts accrued on Board 11 during test
    matrix[N_id, N_food] counts_b11;
    for (j in 1:N_id) {
      counts_b11[j, ] = rep_row_vector(0.0, N_food); // Start both food-items at 0 
    }
    
    // Initialise padded outputs and per-crow trial counters
    array[N_id] int t_ctr;
    for (j in 1:N_id) {
      t_ctr[j] = 0;
      for (f in 1:N_food) {
        for (tm in 1:T_max) {
          a_s[j, f, tm] = 0.0;  // Padding, overwritten if true trial
          a_b[j, f, tm] = 0.0;  // Padding, overwritten if true trial
        }
      }
    }

    // Initialise OH from sampling phase and board food beliefs
    for (j in 1:N_id) {
      alpha_s[j]  =  S_alpha[j];                
      beta_s[j]   = S_beta[j];     
      theta_b11[j] = alpha_s[j] .* beta_s[j]; 
      theta_b12[j] = alpha_s[j] .* beta_s[j]; 
    }

    // Initialise attractions from preference + sampling values
    matrix[N_id, N_food] ALL_F_run;
    for (j in 1:N_id) {
      for (f in 1:N_food) {        
        ALL_F_run[j, f] = inv_logit(ALL_F_value[j, f]); // Raw -> outcome scale
      }
    }

    // Sample-by-sample updates in test
    for (s in 1:N_test) {
      // Loop index object
      int j = id[s];

      // Dirichlet-consistent (i.e., normalised) probabilities matching sampling phase logic
      // Could try softmax and see if changes result
      vector[N_food] P_F_b11 = theta_b11[j] / sum(theta_b11[j]);
      vector[N_food] P_F_b12 = theta_b12[j] / sum(theta_b12[j]);

      // Attractions snapshot
      // No longer sending to global; if want, change back
      vector[N_food] A_F;
      for (f in 1:N_food) { 
        A_F[f] = ALL_F_run[j,f];
      } 

      // Board utilities - prob of food * how much each food valued
      real U_b11 = dot_product(P_F_b11, A_F);
      real U_b12 = dot_product(P_F_b12, A_F);

      // Sensitivity and temporal bias applied to switch only i.e., biasing toward/away from switching
      real L = exp(lambda + v_ID[j,3]); // Raw -> outcome, including offset
      vector[2] util;
      util[1] = L * U_b11;
      util[2] = L * U_b12;
      util[2] += trial_effect * trial[s]; // Asymmetric bias

      // Stay/switch probabilities
      Prob_S[s] = softmax(util);
      
      // Board 11 food counts update
      counts_b11[j, y_test[s]] += 1.0;

      // Overhypothesis updates
      real rho_a = inv_logit(rho_alpha + v_ID[j,1]); // Raw -> outcome, including offset
      real rho_b = inv_logit(rho_beta  + v_ID[j,2]); // Raw -> outcome, including offset
      
      vector[N_food] obs = rep_vector(0.0, N_food); // What food observed eaten?
      obs[y_test[s]] = 1.0;                         // Assign it to correct column i.e., obs[0, 1]
      
      vector[N_food] alpha_new = alpha_s[j] + obs;                 // Pseudo-prior update
      alpha_s[j] = (1.0 - rho_a) * alpha_s[j] + rho_a * alpha_new; // Scaled pseudo-prior update, based on how readily the crow replaces 'old' with 'new' (rho_a)

      vector[N_food] beta_new = alpha_s[j] / sum(alpha_s[j]);   // Pseudo-prior update: we convert alpha_s (pseudo-prior = belief strength)
                                                                // into proportions so that beta_s represents the relative composition of foods.
                                                                // Example: if alpha_s = [4, 6], then beta_new = [0.4, 0.6],
                                                                // meaning the crow expects 40% high-quality and 60% low-quality items.
                                                                // This matches the logic of beta, which always ranges from 0 to 1 and 
                                                                // represents simplex values of belief.  
      beta_s[j] = (1.0 - rho_b) * beta_s[j] + rho_b * beta_new; // Scaled pseudo-prior update, based on how readily the crow replaces 'old' with 'new' (rho_b)

      // Recompute board pseudo-counts
      theta_b11[j] = alpha_s[j] .* beta_s[j] + (counts_b11[j,]'); // OH + Board 11 evidence; crow carries forward level 2 and level 1 knowledge
      theta_b12[j] = alpha_s[j] .* beta_s[j];                     // OH only; crow can only carry forward level 2 knowledge because there's no Board 12 evidence

      // Update attractions for chosen food
      real P = inv_logit(phi + v_ID[j,4]); // Raw -> outcome, including offset
      for (f in 1:N_food) {         
        if (y_test[s] == f) { // Should be f == 2 until switch
          if (y_test_ate[s] + 1 == f) { // If crow ate the low-quality food i.e., 1 + 1 == 2
          // What was sampled updates
            ALL_F_run[j, f] = (1.0 - P) * ALL_F_run[j,f] + P * .25; // This is one option where the low-value item is assigned a low-value score between 0 and 1 
         // ALL_F_run[j, f] = (1.0 - P) * ALL_F_run[j,f] + P * mean(inv_logit(ALL_F_value[, f])); // This is another option where the low-value item is assigned the mean low-quality food attraction (range 0 to 1) across crows from pref + sampling
         // There are other options... 
         
         // Also possible - and typical - to update what wasn't sampled, but it's also okay not to
         // ALL_F_run[j, f - 1] = (1.0 - P) * ALL_F_run[j, f - 1]; // Whatever isn't sampled always gets no reward; f - 1 indexing assumes 2 always coded as low-quality item, which is the only one able to be sampled in test
          
          } else { // If crow did not eat the low-quality food i.e., 0 + 1 == 1 which != 2
              ALL_F_run[j, f] = (1.0 - P) * ALL_F_run[j,f]; // No reward
              
              // And again, need to decide whether not-eaten food attraction updates as well; if it does:
              // ALL_F_run[j, f - 1] = (1.0 - P) * ALL_F_run[j, f - 1]; // Whatever isn't sampled always gets no reward; f - 1 indexing assumes 2 always coded as low-quality item, which is the only one able to be sampled in test
          }
        }
      }
      
      // Save snapshots to padded arrays 
      t_ctr[j] += 1;
      int tloc = t_ctr[j];
      for (f in 1:N_food) {
        a_s[j, f, tloc] = alpha_s[j][f];
        a_b[j, f, tloc] = beta_s[j][f];
      }

      
    }
  }
}

model{
  
  // Preference test priors
  PT_alpha ~ normal(0, 1);                
  PT_weight ~ exponential(1);           
  to_vector(PT_F_value) ~ normal(0, 1); // Vectorise because matrix

  // Sampling priors
  S_gamma ~ exponential(1); // Over-over hypothesis
  for (j in 1:N_id) {
    S_alpha[j] ~ exponential(S_gamma[j]);             // Overhypothesis about food item variability
    S_beta[j]  ~ dirichlet(rep_vector(1.0, N_food));  // Overhypothesis about food item distribution
    for (board in 1:N_board) {
      S_theta[j, board] ~ dirichlet(S_alpha[j] .* S_beta[j]);
    }
  }
  
  // Random effects priors
  to_vector(z_ID) ~ normal(0, 1);
  sigma_ID ~ exponential(1);
  Rho_ID ~ lkj_corr_cholesky(4);
  
  //Sampling likelihood  
  for(j in 1:N_id){
    for(board in 1:N_board){
      y_sampling[j, board] ~ multinomial(S_theta[j, board]);
    }
  }
  
  //Test priors...
  
  //Attraction priors
  ALL_alpha ~ normal(0, 1);                
  ALL_weight ~ exponential(1);              
  to_vector(ALL_F_value) ~ normal(0, 1);
  
  //Weighting parameters
  //Higher values expected for earlier switchers = implies rapid adaptation to test-phase data
  //Lower values expected for later switchers = implies stronger reliance on prior knowledge 
  rho_alpha ~ normal(0, 1); 
  rho_beta ~ normal(0, 1);
  phi ~ normal(0, 1); 

  //Scaling parameter
  //Higher values expected for earlier switchers = implies increased sensitivity to perceived probabilistic diff's in food-item probabilities per unmeasured factors
  //Lower values expected for later switchers = less fussed about perceived probabilistic diff's in food-item probabilities per unmeasured factors
  //Such unmeasured factors might be (but are not limited to):
  //1) Cost of switching e.g., if viewed as high, switch less likely despite accumulating 'bad' food evidence (lambda would be low)
  //2) Frustration with repeated 'bad' food samples e.g., if easily frustrated, more likely to readily switch (lambda would be high)
  lambda ~ normal(0, 1); 
  
  //Time parameter
  //Positive values expected for earlier switchers = implies rapid keeness to switch
  //Negative values expected for later switchers = implies growing resistance to switching over time
  //Such keeness or resistance could arise e.g., via 'personality' differences or familiarity (with board) effects, respectively
  trial_effect ~ normal(0, 1);

  //NOTE: having both the time & scaling parameters allows for the disentangling of effects
  //e.g., could have low sensitivity to diff's in food-item probabilities (L) but strong positive temporal effect (early switch)
  
  //Determine...
  
  //Food attractions in pref test and across pref test plus sampling
  for(j in 1:N_id){
    for(f in 1:N_food){
      PT_total_eat[j,f] ~ poisson_log(PT_alpha + PT_log_dur[j] + PT_weight * inv_logit(PT_F_value[j, f]));
      ALL_total_eat[j,f] ~ poisson_log(ALL_alpha + ALL_log_dur[j] + ALL_weight * inv_logit(ALL_F_value[j, f]));
    }
  }
  
  //Switch likelihood
  for(s in 1:N_test){
    Switch[s] ~ categorical(Prob_S[s]);  
  }
}

