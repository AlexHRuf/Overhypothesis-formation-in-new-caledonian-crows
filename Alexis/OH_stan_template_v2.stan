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
}
parameters{
  //Sampling parameters
  real<lower=0> gamma[N_id]; //Over-over hypothesis
  vector<lower=0>[N_colors] alpha[N_id]; //Overhypthosis about food item variability
  simplex[N_colors] beta[N_id]; //Overhypothesis about food item distribution
  array[N_id] simplex[N_food] theta[N_board]; //Local knowledge about food item variability & distribution
  
  //Test parameters
  real<lower=0, upper=1> weight_alpha[N_id]; //Weight for updating alpha dynamically in-test
  real<lower=0, upper=1> weight_beta[N_id]; //Weight for updating beta dynamically in-test
}
model{
  //Sampling priors
  for(j in 1:N_id){
    gamma[j] ~ exponential(1);
    alpha[j] ~ exponential(gamma);
    beta[j] ~ dirichlet(rep_vector(1.0, N_food));
    for(board in 1:N_board){
      theta[j][board] ~ dirichlet(alpha[j] .* beta[j]);
    }
  }
  //Sampling likelihood
  for(j in 1:N_id){
    for(board in 1:N_board){
      y_sampling[j][board, :] ~ multinomial(theta[j][board]); 
    }
  }
  
  //Test priors
  weight_alpha ~ normal(0, 1);
  weight_beta ~ normal(0, 1);
  
  //Define & initialise objects for sequential updating for test (board 11)
  matrix[N_id, 2] food_counts; //Running food count per item per bird
  vector<lower=0>[N_food] alpha_s[N_id]; //Step-specific alpha for each individual
  simplex[N_food] beta_s[N_id]; //Step-specific beta for each individual

  for(j in 1:N_id){
    food_counts[j, 1:2] = rep_vector(0.0, N_food)';
    alpha_s[j] = alpha[j];
    beta_s[j] = beta[j];
  }

  //Loop through sample choices per individual
  for(s in 1:N_test){
    //Prob of food item
    vector[N_colors] P_F; //Probabilities of each food item before sampling
    real w_alpha; 
    real w_beta; 
    vector[N_colors] new_alpha;
    //Update posterior counts based on observed sample
    if (s > 1){
      food_counts[id[s], y_test[s - 1]] += 1.0; //Increment count for observed sample
    }
    //Compute new alpha based on posterior counts
    new_alpha = alpha[id[s]] + food_counts[id[s], :]; 
    //Weighted update of alpha
    w_alpha = inv_logit(weight_alpha[id[s]]); //THO MAYBE DO COV MATRIX FOR ID VAR
    alpha_s[id[s]] = w_alpha * new_alpha + (1 - w_alpha) * alpha_s[id[s]];
    //Update beta dynamically based on new alpha
    w_beta = inv_logit(weight_beta[id[s]]); //THO MAYBE DO COV MATRIX FOR ID VAR
    beta_s[id[s]] = w_beta * (alpha_s[id[s]] / sum(alpha_s[id[s]])) + (1 - w_beta) * beta_s[id[s]];
    //Compute posterior probabilities for the next potential food sample
    P_F = softmax(food_counts[id[s], :] + alpha_s[id[s]] .* beta_s[id[s]]);
  }
}
