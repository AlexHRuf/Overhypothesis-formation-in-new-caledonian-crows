data{

  //Sampling variables
  int<lower=1> N_id; //Total number of individuals
  int<lower=1> N_board; //Number of boards in the sampling phase (10)
  int<lower=1> N_food; //Number of food items (2)
  array[N_id] matrix[N_board, N_food] y_sampling; //Total counts of each food item per board per bird
  
  //Test variables
  int<lower=1> N_test; //Number of observed samples from board 11
  int<lower=1> y_test[N_test]; //Observed samples from bag 11 (1 = high quality item; 2 = low quality item)
  int<lower=1> id[N_test]; //Unique identification
  int<lower=1> trial_num[N_test]; // Trial number (e.g., 1, 2, 3, ...)
  int<lower=1> Switch[N_test]; //1 = no; 2 = yes
  
}
parameters{

  //Sampling parameters
  real<lower=0> gamma[N_id]; //Over-over hypothesis
  vector<lower=0>[N_food] alpha[N_id]; //Overhypthosis about food item variability
  simplex[N_food] beta[N_id]; //Overhypothesis about food item distribution
  array[N_id] simplex[N_food] theta[N_board]; //Local knowledge about food item variability & distribution
  
  //Test parameters
  real weight_alpha[N_id]; //Weight for updating alpha dynamically in-test
  real weight_beta[N_id]; //Weight for updating beta dynamically in-test
  real lambda[N_id]; //Sensitivity about P_F differences
  real trial_effect; // Effect of trial number on switch likelihood
}
model{

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
  
  //Test priors  
  //Weighting parameters
  //Higher values expected for earlier switchers = implies rapid adaptation to test-phase data
  //Lower values expected for later switchers = implies stronger reliance on prior knowledge 
  weight_alpha ~ normal(0, 1); 
  weight_beta ~ normal(0, 1);

  //Scaling parameter
  //Higher values expected for earlier switchers = implies increased sensitivity to perceived probabilistic diff's in food-item probabilities per unmeasured factors
  //Lower values expected for later switchers = less fussed about perceived probabilistic diff's in food-item probabilities per unmeasured factors
  //Such unmeasured factors might be (but are not limited to):
  //1) Cost of switching e.g., if viewed as high, switch less likely despite accumulating 'bad' food evidence (lambda would be low)
  //2) Frustration with repeated 'bad' food samples e.g., if easily frustrated, more likely to readily switch (lambda would be high)
  Lambda ~ normal(0, 1); 

  //Time parameter
  //Positive values expected for earlier switchers = implies rapid keeness to switch
  //Negative values expected for later switchers = implies growing resistance to switching over time
  trial_effect ~ normal(0, 1);

  //Note: having both the time & scaling parameters allows for the disentangling of effects
  //e.g., could have low sensitivity to diff's in food-item probabilities (L) but strong positive temporal effect (early switch)
  
  //Define objects for sequential updating for test (board 11)
  matrix[N_id, 2] food_counts; //Running food count per item per bird
  vector<lower=0>[N_food] alpha_s[N_id]; //Step-specific alpha for each individual
  simplex[N_food] beta_s[N_id]; //Step-specific beta for each individual

  //Initialise these objects
  for(j in 1:N_id){
    food_counts[j, :] = rep_vector(0.0, N_food)'; //No food yet sampled
    alpha_s[j] = alpha[j]; //Initialise to estimates over sampling phase
    beta_s[j] = beta[j]; //Initialise to estimates over sampling phase
  }

  //Loop through sample choices per individual
  for(s in 1:N_test){
    //Local loop objects
    vector[N_food] P_F; //Probabilities of food item
    vector[N_food] P_2; //Probability of switch
    real w_alpha; 
    real w_beta; 
    real L; //Lambda
    vector[N_food] new_alpha;
    vector[N_food] new_beta;
    real trial_influence;
    
    //Update food counts based on observed sample
    if (s > 1){
      food_counts[id[s], y_test[s - 1]] += 1.0; //Increment count for observed sample
    }
    
    //Compute new alpha based on posterior counts
    new_alpha = alpha[id[s]] + food_counts[id[s], :]; 
    
    //Weighted dynamic update of alpha
    w_alpha = inv_logit(weight_alpha[id[s]]); //THO MAYBE DO COV MATRIX FOR ID VAR
    alpha_s[id[s]] = w_alpha * new_alpha + (1 - w_alpha) * alpha_s[id[s]];
    
    //Update beta dynamically based on new alpha
    new_beta = alpha_s[id[s]] / sum(alpha_s[id[s]]); //Normalise alpha_s (vectorised -> simplex)
    w_beta = inv_logit(weight_beta[id[s]]); //THO MAYBE DO COV MATRIX FOR ID VAR
    beta_s[id[s]] = w_beta * new_beta + (1 - w_beta) * beta_s[id[s]];
    
    //Compute food probabilities 
    P_F = softmax(food_counts[id[s], :] + alpha_s[id[s]] .* beta_s[id[s]]);

    //Initialise sensitivity scaler
    L = exp(Lambda[id[s]]);

    //Compute switch probabilities
    trial_influence = trial_effect * trial_num[s]; // Scale trial numbers
    P_S = softmax(L * P_F + trial_influence);

    //Switch likelihood
    Switch[s] ~ categorical(P_S);  
    
  }
}
