data{

  //General variables
  int<lower=1> N_test; //Number of overall observations in test phase
  int<lower=1> N_food; //Number of food items (2)
  int<lower=1> N_id;   //Total number of individuals
  
  //Preference variables
  int PT_total_eat[N_id, N_food]; //Total pieces of food item 1 (high) and 2 (low) eaten
  real PT_log_dur[N_id, N_food];  //Total duration of preference test in seconds (log transformed)

  //Sampling variables
  int<lower=1> N_board;                           //Number of boards in the sampling phase (10)
  array[N_id] matrix[N_board, N_food] y_sampling; //Total counts of each food item per board per bird
  
  //Test variables
  int ALL_total_eat[N_id, N_food]; //Total pieces of food item 1 and 2 eaten across pref test & sampling per bird
  real ALL_log_dur[N_id];          //Total duration (log transformed) spent across pref test & sampling per bird
  int<lower=1> y_test[N_test];     //Observed samples from board 11 (1 = high quality item; 2 = low quality item)
  int<lower=1> y_test_ate[N_test]; //Observed eating (or not) of sample from board 11
  int<lower=1> id[N_test];         //Unique identification
  int<lower=1> trial[N_test];      //Trial number (e.g., 1, 2, 3, ...)
  int<lower=1> Switch[N_test];     //1 = no; 2 = yes
  
}

parameters{
  
  //Preference test latent parameters
  real PT_alpha;                   //Intercept i.e., rate of random choice 
  real<lower=0> PT_weight;         //Weight for influence of PT_F_value on total choices; some birds might be more affected than others
  matrix[N_id, N_food] PT_F_value; //High and Low quality individual-level 'value' estimates for all food eaten in pref test

  //Sampling latent parameters parameters
  real<lower=0> gamma[N_id];                  //Over-overhypothesis
  vector<lower=0>[N_food] alpha[N_id];        //Overhypothesis about food item variability
  simplex[N_food] beta[N_id];                 //Overhypothesis about food item distribution
  array[N_id] simplex[N_food] theta[N_board]; //Local knowledge about food item variability & distribution
  
  //Test latent parameters
  real ALL_alpha;                   //Intercept i.e., rate of random choice 
  real<lower=0> ALL_weight;         //Weight for influence of ALL_F_value on total choices; some birds might be more affected than others
  matrix[N_id, N_food] ALL_F_value; //High and Low quality individual-level 'value' estimates for all food eaten across pref test & sampling
  real weight_alpha;                //Weight for updating alpha dynamically in-test
  real weight_beta;                 //Weight for updating beta dynamically in-test
  real lambda;                      //Sensitivity about P_F differences in boards
  real phi;                         //Learning rate about food attractions
  real trial_effect;                //Effect of trial number on switch likelihood
  
  //Individual varying effects parameters
  matrix[4, N_id] z_ID;           
  vector<lower = 0>[4] sigma_ID;  
  cholesky_factor_corr[4] Rho_ID; 
  
}

transformed parameters{
  
  //Global scope i.e., saved in output
  matrix[N_id, 4] v_ID;  
  simplex[N_food] Prob_S[N_test]; //Prob switch for posterior prediction checks
  matrix[N_test, N_food] a_s;     //Step-specific alpha 
  matrix[N_test, N_food] a_b;     //Step-specific beta
  matrix[N_test, N_food] A;       //Step-specific attractions
  
  v_ID = (diag_pre_multiply(sigma_ID, Rho_ID) * z_ID)'; 
  
  //NOTE: It is possible to pull/graph above final four global scope params by bird, given that we know max test trials per bird from raw data
  //e.g., if bird 1 did 10 test trials & bird 2 did 15 test trials
  //a_s[,1:10,] will be for bird 1 & a_s[,11:25,] will be for bird 2, when accessing posterior samples
  
  { //Local scope i.e., not saved in output
  
    //Define objects for sequential updating for test (board 11)
    vector<lower=0>[N_food] alpha_s[N_id]; //Step-specific alpha for each individual
    simplex[N_food] beta_s[N_id];          //Step-specific beta for each individual
    matrix[N_id, 2] food_counts;           //Running food count per item per bird
    matrix[N_id, 2] ALL_F;                 //Running food attractions per item per bird
    
    //Initialise...
    
    //Attraction to food across pref test & sampling 
    for(j in 1:N_id){
      for(f in 1:N_food){
        ALL_F[j, f] = inv_logit(ALL_F_value[j, f]); 
      }
    }
    
    //Running food counts (to zero) & dynamic alpha and beta (to sampling phase estiamtes)
    for(j in 1:N_id){
      food_counts[j, :] = rep_vector(0.0, N_food)'; 
      alpha_s[j] = log(alpha[j]); 
      beta_s[j] = beta[j]; 
    }
  
    //Loop through sample choices per individual
    for(s in 1:N_test){
      
      //Local loop objects
      vector[N_food] P_F;       //Probabilities of high and low food items
      vector[N_food] A_F;       //Attractions to high and low food items
      vector[N_food] P_S;       //Probabilities of stay or switch
      vector[N_food] new_alpha; //Unweighted dynamic alpha
      vector[N_food] new_beta;  //Unweighted dynamic beta
      real w_alpha;             //Alpha weight
      real w_beta;              //Beta weight
      real L;                   //Lambda - sensitivity to probabalistic diff's in board food-payoffs
      real P;                   //Phi - food attraction updating rate
      real trial_influence;     //Temporal trial effect
      
      //Compute new alpha based on posterior counts
      new_alpha = alpha[id[s]] + food_counts[id[s], :]; 
      
      //Weighted dynamic update of alpha
      w_alpha = inv_logit(weight_alpha + v_ID[id[s], 1]); 
      alpha_s[id[s]] = w_alpha * new_alpha + (1 - w_alpha) * alpha_s[id[s]];
      a_s[s, :] = alpha_s[id[s]]'; //To global
      
      //Update beta dynamically based on new alpha & give appropriate weighting
      new_beta = alpha_s[id[s]] / sum(alpha_s[id[s]]); //Normalise alpha_s (vectorised -> simplex)
      w_beta = inv_logit(weight_beta + v_ID[id[s], 2]); 
      beta_s[id[s]] = w_beta * new_beta + (1 - w_beta) * beta_s[id[s]];
      a_b[s, :] = beta_s[id[s]]'; //To global

      //Compute food probabilities for boards 11 & 12
      P_F_b11 = softmax(food_counts[id[s], :] + alpha_s[id[s]] .* beta_s[id[s]]); //Responds to local sample experience
      P_F_b12 = softmax(alpha_s[id[s]] .* beta_s[id[s]]);                         //Reflects only overhypothesis
  
      //Initialise sensitivity scaler, food attractions & teporal influence
      L = exp(lambda + v_ID[id[s], 3]); 
      A_F = ALL_F[id[s], :]; //Matrix -> vector so easier calcs in P_S
      A[s, :] = ALL_F[id[s]]; //To global
      trial_influence = trial_effect * trial[s]; 

      //Compute switch probabilities
      P_S = softmax(L * ( (P_F_b11 * A_F) / (P_F_b12 * A_F) ) + trial_influence);
      Prob_S[s] = P_S; //Send to global
      
      //Update food counts after observed sample
      food_counts[id[s], y_test[s]] += 1.0; 
      
      //Update attraction scores after observed sample
      P = inv_logit(phi + v_ID[id[s], 4]); //THO MAYBE DO COV MATRIX FOR ID VAR
      for(f in 1:N_food){
        if(y_test[s] == f){
          ALL_F[id[s], f] = (1 - P) * ALL_F[id[s], f] + P * y_test_ate[s]; //Update only for sampled item
        } 
      }
      
      //NOTE: have to think about whether y_test_ate (acting as the pay variable) as 0 (no) vs 1 (yes) is best parameterisation
      //Because if eats the low value item, attraction will increase for that item
      //One idea could be to have low quality item reward something (e.g., 0.25 or 0.5) if eaten, and 0 if not eaten
      //Alternatively, this might be where we could have y_test_ate pay according to individual value estimates from pref test per food type i.e., inv_logit(PT_F_value[j, f])
      //It all depends on how often the low value item is generally eaten from board 11 - I don't know this
      //WHAT DO YOU THINK ALEX?
      
    } //End s loop
  } //End local scope
  
}

model{
  
  //Preference test priors
  PT_alpha ~ normal(0, 1);                
  PT_weight ~ exponential(1);              
  to_vector(PT_F_value) ~ normal(0, 1);

  //Sampling priors
  for(j in 1:N_id){
    gamma[j] ~ exponential(1); //Over-over hypothesis
    alpha[j] ~ exponential(gamma); //Overhypothesis about food item variability
    beta[j] ~ dirichlet(rep_vector(1.0, N_food)); //Overhypothesis about food item distribution
    for(board in 1:N_board){
      theta[j][board] ~ dirichlet(alpha[j] .* beta[j]); //Local knowledge about food item variability & distribution
    }
  }
  
  //Sampling likelihood  
  for(j in 1:N_id){
    for(board in 1:N_board){
      y_sampling[j][board, :] ~ multinomial(theta[j][board]); 
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
  weight_alpha ~ normal(0, 1); 
  weight_beta ~ normal(0, 1);
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


