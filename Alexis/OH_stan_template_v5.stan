data{
  // General indexing 
  int<lower=1> N_test;        // Total test-phase observations 
  int<lower=1> N_food;        // Total food items (2)
  int<lower=1> N_id;          // Total individuals (16?)
  int<lower=1> N_board;       // Total boards in sampling phase (10)
  int<lower=1> T_id[N_id];    // Trials per crow in test (can be 0 for some)
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
  vector<lower=0>[N_id] gamma_s;          // Over-overhypothesis
  vector<lower=0>[N_food] alpha_s[N_id];  // Overhypothesis about food item variability
  simplex[N_food] beta_s[N_id];           // Overhypothesis about food item distribution
  simplex[N_food] theta_s[N_id, N_board]; // Local (per board) knowledge about food item variability & distribution

  // Pref + sampling attractions 
  real ALL_alpha;                   // Intercept i.e., rate of random choice 
  real<lower=0> ALL_weight;         // Weight for food-item attraction values i.e., allow attractions to differentially direct eating, because can't assume 'like' or 'dislike' has same effect on food-choice across crows
  matrix[N_id, N_food] ALL_F_value; // High and Low quality individual-level attraction value estimates for all food eaten 

  // Test phase target parameters on raw scale
  real lambda_raw;               // Sensitivity rate (weight) to differences in board utility (i.e., prob food informed by OH * attraction to food)
  real phi_raw;                  // Food attraction values updating rate (scaler)
  vector[N_food] alpha_last_raw; // Endpoint alpha
  real beta_last_raw;            // Endpoint beta (of low-quality food distribution prior to transformation - see below)
  
  // Test phase monotonic effects
  matrix[N_id, T_max] delta_alpha_raw;  
  matrix[N_id, T_max] delta_beta_raw;   
  
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
  array[N_id, N_food, T_max] real a_t; // Padded snapshot OH alpha_t (i.e., zeros in unused slots); facilitates array indexing b/c otherwise would be ragged structure
  array[N_id, N_food, T_max] real b_t; // Padded snapshot OH beta_t  (i.e., zeros in unused slots); facilitates array indexing b/c otherwise would be ragged structure

  { // Local scope - nothing below here saved
  
    // Dynamic variable-updates per crow in-test
    array[N_id] vector[N_food] alpha_t;   // OH variability based on monotonic effects structure
    array[N_id] vector[N_food] beta_t;    // OH distribution based on monotonic effects structure
    array[N_id] vector[N_food] theta_b11; // Board 11 belief = OH * local knowledge
    array[N_id] vector[N_food] theta_b12; // Board 12 belief = OH only
    
    // Initialise food counts accrued on Board 11 during test
    matrix[N_id, N_food] counts_b11;
    for (j in 1:N_id) {
      counts_b11[j, ] = rep_row_vector(0.0, N_food); // Start both food-items at 0 
    }
    
    // Initialise per-crow trial counters and padded outputs
    array[N_id] int t_ctr;
    for (j in 1:N_id) {
      t_ctr[j] = 0;
      for (f in 1:N_food) {
        for (tm in 1:T_max) {
          a_t[j, f, tm] = 0.0; // Padding, overwritten if true trial
          b_t[j, f, tm] = 0.0; // Padding, overwritten if true trial
        }
      }
    }

    // Initialise first parameters, dynamic parameters & board food-probabilities
    array[N_id] vector[N_food] a_first;
    array[N_id] vector[N_food] b_first;
    for (j in 1:N_id) {
      a_first[j] = alpha_s[j];                // Alpha first always OH about food variability prior to test feedback
      b_first[j] = beta_s[j];                 // Beta first always OH about food distribution prior to test feedback
      alpha_t[j]  = a_first[j];               // Dynamic alpha initialised to alpha first prior to to test feedback
      beta_t[j]   = b_first[j];               // Dynamic beta initialised to beta first prior to test feedback
      theta_b11[j] = alpha_t[j] .* beta_t[j]; // OH across Boards 1 - 10
      theta_b12[j] = alpha_t[j] .* beta_t[j]; // OH across Boards 1 - 10
    }

    // Initialise attractions from preference + sampling values
    matrix[N_id, N_food] ALL_F_run;
    for (j in 1:N_id) {
      for (f in 1:N_food) {        
        ALL_F_run[j, f] = inv_logit(ALL_F_value[j, f]); // Raw -> outcome scale
      }
    }
    
    // Softmax-prefix deltas (make per-bird, per-trial weights that sum to 1) 
    // Goal: turn unconstrained raw deltas (N_id × T_max) into, for each bird j,
    // a length-T_id[j] vector of nonnegative weights that sum to 1.
    // These weights allocate *how much* of the total change happens at each trial.
    // We store them in delta_param_w; columns beyond T_id[j] stay 0 as padding.

    matrix[N_id, T_max] delta_alpha_w = rep_matrix(0, N_id, T_max);  // Start with all zeros
    matrix[N_id, T_max] delta_beta_w = rep_matrix(0, N_id, T_max);   // Unused tail will remain 0

    for (j in 1:N_id) {
      // 1) Take the j-th row of the raw deltas, convert to a column vector,
      // 2) Keep only the first T_id[j] entries (this bird's actual trials),
      // 3) Map through softmax to get weights greater than or equal to 0 that sum to 1 across those trials.
      //
      // Notes:
      // - to_vector(...) converts the row to a column vector (softmax expects a vector).
      // - head(x, k) slices the first k entries; this discards any extra T_max - T_id[j] tail.
      // - softmax(...) exponentiates and normalises so the k entries sum to 1.
      vector[T_id[j]] wa = softmax(head(to_vector(delta_alpha_raw[j]), T_id[j]));
      vector[T_id[j]] wb = softmax(head(to_vector(delta_beta_raw[j]), T_id[j]));
    
      // 4) Write those weights back into columns 1:T_id[j] of the j-th row.
      //    We transpose (') to turn column vector into row vector expected by delta matrixes
      //    Columns T_id[j]+1 : T_max remain 0 (padding for rectangular arrays).
      delta_alpha_w[j, 1:T_id[j]] = wa';
      delta_beta_w [j, 1:T_id[j]] = wb';
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

      // Sensitivity to board utilities 
      real L = exp(lambda_raw + v_ID[j,1]); // Raw -> outcome, including offset
      
      // Apply sensitivity to board utilities, and include any temporal bias (applied to switch only i.e., biasing toward/away from switching)
      vector[2] util;
      util[1] = L * U_b11;
      util[2] = L * U_b12;
      util[2] += trial_effect * trial[s]; // Asymmetric bias

      // Stay/switch probabilities
      Prob_S[s] = softmax(util);
      
      // Board 11 food counts update
      counts_b11[j, y_test[s]] += 1.0;

      // Overhypothesis updates
      t_ctr[j] += 1; // Update per bird trial counter to correctly index into deltas
      
      // Current deltas 
      real cumA = sum(delta_alpha_w[j, 1:t_ctr[j]); 
      real cumB = sum(delta_beta_w[j, 1:t_ctr[j]); 
      
      // Last params
      vector[2] a_last = exp(alpha_last_raw + rep_vector(v_ID[j, 2], 2)); // Raw -> outcome, including offset; repped 2 x to match vector structure of a_last
      
      vector[2] b_last;
      real b_j_last_raw = beta_last_raw + v_ID[j, 3]; // Draw from raw & add raw offset
      b_last[2] = inv_logit(b_j_last_raw);            // Raw -> outcome, arbitarily set to refer to distribution on low-quality food item
      b_last[1] = 1 - b_last[2];                      // Subtracting that from 1 makes a simplex, consistent with beta logic in sampling 

      // Current alpha and beta
      alpha_t[j] = a_first[j] + (a_last - a_first[j]) * cumA;
      beta_t[j] = b_first[j] + (b_last - b_first[j]) * cumB; 

      // Recompute board pseudo-counts mirror sampling logic
      vector[N_food] OH = alpha_t[j] .* beta_t[j];
      theta_b11[j] = OH + (counts_b11[j, ]') ;
      theta_b12[j] = OH;

      // Update attractions for chosen food
      real P = inv_logit(phi_raw + v_ID[j,4]); // Raw -> outcome, including offset
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
      for (f in 1:N_food) {
        a_t[j, f, t_ctr[j]] = alpha_t[j][f];
        b_t[j, f, t_ctr[j]] = beta_t[j][f];
      }

      
    }
  }
}

model{
  
  // Preference test priors
  PT_alpha ~ normal(0, 1);                
  PT_weight ~ exponential(1);           
  to_vector(PT_F_value) ~ normal(0, 1); // Vectorise because matrix

  // Sampling phase priors
  gamma_s ~ exponential(1); // Over-over hypothesis
  for (j in 1:N_id) {
    alpha_s[j] ~ exponential(gamma_s[j]);             // Overhypothesis about food item variability
    beta_s[j]  ~ dirichlet(rep_vector(1.0, N_food));  // Overhypothesis about food item distribution
    for (board in 1:N_board) {
      theta_s[j, board] ~ dirichlet(alpha_s[j] .* beta_s[j]);
    }
  }
  
  //Sampling likelihood  
  for(j in 1:N_id){
    for(board in 1:N_board){
      y_sampling[j, board] ~ multinomial(theta_s[j, board]);
    }
  }
  
  //Test priors...
  
  // Random effects priors
  to_vector(z_ID) ~ normal(0, 1);
  sigma_ID ~ exponential(1);
  Rho_ID ~ lkj_corr_cholesky(4);
  
  //Attractions
  ALL_alpha ~ normal(0, 1);                
  ALL_weight ~ exponential(1);              
  to_vector(ALL_F_value) ~ normal(0, 1);
  
  // Monotonic effects 
  alpha_last_raw ~ normal(0, 1);                  
  beta_last_raw  ~ normal(0, 1);
  to_vector(delta_alpha_raw) ~ normal(0, 1);          
  to_vector(delta_beta_raw)  ~ normal(0, 1);
  
  //Scaling parameter
  //Higher values expected for earlier switchers = implies rapid adaptation to test-phase data
  //Lower values expected for later switchers = implies stronger reliance on prior knowledge 
  phi_raw ~ normal(0, 1); 

  //Weighting parameter
  //Higher values expected for earlier switchers = implies increased sensitivity to perceived probabilistic diff's in food-item probabilities per unmeasured factors
  //Lower values expected for later switchers = less fussed about perceived probabilistic diff's in food-item probabilities per unmeasured factors
  //Such unmeasured factors might be (but are not limited to):
  //1) Cost of switching e.g., if viewed as high, switch less likely despite accumulating 'bad' food evidence (lambda would be low)
  //2) Frustration with repeated 'bad' food samples e.g., if easily frustrated, more likely to readily switch (lambda would be high)
  lambda_raw ~ normal(0, 1); 
  
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
