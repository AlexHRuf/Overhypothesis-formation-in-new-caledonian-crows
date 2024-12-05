data{ 

  int N;           //Total number of observations       
  int N_id;        //Total number of birds      
  int N_food;      //Total number of food items (2)          
  int N_sess;      //Total number of sessions       
  int T_test;      //Total number of test trials across individuals        
  
  int total_eat[N_id, N_food]; //Total pieces of food item 1 and 2 eaten
  real log_dur[N_id];          //Total duration of preference test in seconds (log transformed)
  //Count time to for preferences for each bird. Sum entire timeframe of exposure.
  
  int id[N];        //Bird-specific individual identificaiton number (1 .. N)        
  int sess[N];      //Session identification number (1 .. N_sess)           
  int trial[N];     //Within session running trial account      
  int cum_trial[N]; //Running trial count that resets between phases      
  int food_sample[N];    //Food item sampled (1 = high, 2 = low)
  int ate_sample[N];     //Wheter the bird ate the food item (0 = no, 1 = yes)  
  int test[N];      //Whether the test session (0 = no, 1 = yes)   
  int choice[T_test];    //Whether stay or switch (1 = stay, 2 = switch) - note: cannot 0/1 dummy code dependent variabels in stan    

}
  
parameters{
    
  //Parameters to model food preference (or value)
  matrix[N_id, N_food] V_p; //Baseline value matrix for each food type and ID
  real<lower=0> weight_p; //Weight parameter: Controls influence of V_p on the total food consumption (total_eat).
  real alpha_p; //Global interecept to account for individual differences in food consumption
  //Note: Change this to account for the actual structure of food preference testing (Trial-by-trial choices).
  
  // Parameters to model the prior probabilities of sampling food types.
  simplex[2] raw_beta[N_id, N_sess];    // Session/ID specific hyperprior for raw_alpha. Represents global believ about prob of food types.
  simplex[2] raw_alpha[N_id, N_sess];   // Session/ID specific prior for the probability of sampling food types (prob_F).
  matrix[N_id, N_sess] logit_phi;       // Session/ID specific learning parameter. Controls how sampling updates prob_F. Used to calculate phi_ind.
  matrix[N_id, N_sess] log_lambda;      // Session/ID specific scaling parameter for the probability of staying. Used to calculate lambda_ind.

  //Parameters to model individual differences (random effects).
  matrix[3, N_id] z_ID; #Later transformed to capture individual-specific effects.      
  vector<lower = 0>[3] sigma_ID; //Standard deviation
  cholesky_factor_corr[3] Rho_ID; //Correlation matrix of parameters
 
} 

transformed parameters{ 
  //global scope (everything here gets posteriors)
  matrix[N_id, 3] v_ID; //Set up matrix for individual parameters      
  //simplex[2] S_Prob[T_test]; //Set up a simplex for stay probabilities in test
  matrix[N, 5] alpha_T; //Set up a matrix for trial by trial alpha updates
  matrix[N, 5] Prob_F;  //Set up a matrix for trial by trial food type probabilities
  //matrix[T_test, 5] Prob_S;  //Set up a matrix for trial by trial (test trials) stay probabilities

  v_ID = (diag_pre_multiply(sigma_ID, Rho_ID) * z_ID)'; //Calculate individual differences using sigma, rho, z. Will later be used to adjust individual phi, lambda, and alpha
  //Partial pooling (check out in rethinking.)
  { //Local scope from here on out (not saved)

    matrix[N_id, N_food] F; //Matrix containing total food eaten per bird
    vector[2] a_T[N_id]; //Set up alpha vector
    matrix[N_id, N_food] Pay; //Set up Pay "for Payoff" matrix 
    matrix[N_id, N_food] A; //Set up A matrix "A for Attraction" A dynamic way to update food preferences during the sampling process? Needs clarification.
   
    for(j in 1:N_id){
    
      A[j, 1:2] = rep_vector(0.1, N_food)'; //Fill A matrix with 0.1 starting values for each food type / bird. Indicates low initial attraction to both food types
      
      for(f in 1:N_food){//Pay is our pi value in EQ 1 and 2.
        Pay[j, f] = inv_logit(V_p[j, f]); //Fill Pay matrix with inv_logit of the value matrix (0-1 scale for pi.)
      }
    
    }

    for(i in 1:N){ //for each observation
      
      //Adjust session-level parameters by individual random effects. 
      real phi_ind = inv_logit(logit_phi[id[i], sess[i]] + v_ID[id[i], 1]); //Individual learning rate phi (for attraction updating). logit_phi[id[i], sess[i]] is the ID/session specific learning rate, + v_ID[id[i], 1] is the individual adjustment. inv_logit() converts the adjusted log odds to a value between 0 and 1
      real lambda_ind = exp(log_lambda[id[i], sess[i]] + v_ID[id[i], 2]); //Individual scaling factor for probability calculations. ID/session specific log_lambda, adjusted by + v_ID[id[i], 2]
      //real phi_ind = 0.5;
      //real lambda_ind = 1.0;

      vector[2] alpha = raw_alpha[id[i], sess[i]]; //Retrieve ID/session specific prior probabilities of sampling each food type (alpha)
      real alpha_unconstrained = log((alpha[1] + 1e-8) / (alpha[2] + 1e-8)); //Calculate the log-odds of the two probabilities. +1e8 to avoid dividsion by 0. Do this so we can add alpha to the offset.
      real alpha_v_ID_sum = alpha_unconstrained + v_ID[id[i], 3]; //Adjust by the ID specific adjustment parameter for alpha v_ID[id[i], 3]
      vector[2] alpha_v_ID_sum_transformed = softmax([alpha_v_ID_sum, 0]'); //Transform adjusted log-odds to a probability vector using softmax. 0 is a reference element (formatting).
      vector[2] alpha_ind = alpha_v_ID_sum_transformed; //Save as alpha_ind. Define as simplex?
      
      vector[5] temp; //Vector to store temporary values
      vector[2] pF; //probability of sampling each food type during the current trial
      vector[2] A_row;
      vector[2] pF_x_A_row;
      vector[2] scaled_pF_x_A_row;
      vector[2] pS;
      vector[2] a_T_update; //Used to update a_T based on new counts.
      
      //point estimte
      if(trial[i] == 1){
        a_T[id[i]] = alpha_ind; //For ID, initialize alpha for trial 1 as a_T
        F[id[i], 1:2] = rep_vector(0, N_food)'; //For ID, initialze counts of each food type at 0 for trial 1.
      }
      
      temp[1] = id[i];
      temp[2] = sess[i];
      temp[3] = trial[i];
      temp[4] = a_T[id[i]][1];
      temp[5] = a_T[id[i]][2];
      
      alpha_T[i, 1:5] = temp'; //Save information for the current trial (id,sess,trial,a_T[1],a_T[2]) in alpha_T
      
      for(f in 1:N_food){ //for each food type (in each trial)
        pF[f] = (F[id[i], f] + a_T[id[i]][f]) / (trial[i] + N_food); //Calculate the probability of sampling each food type f
        
        //F[id[i], f] is the count of food f that ID has sampled up to the current trial
        //+ a_T[id[i]][f] is the prior belief (probability of sampling food type f)
        // devided by (trial[i] + N_food)
      }
      
      pF = pF / sum(pF); //Normalize. Do these have to sum to 1 here?
      
      // After calculating pF, check for invalid values
      for (f in 1:N_food) {
        if (is_nan(pF[f]) || is_inf(pF[f])) {
          reject("pF[", f, "] is invalid at i =", i, ", id[i] =", id[i], ", trial[i] =", trial[i]);
        }
        if (pF[f] < 0 || pF[f] > 1) {
          reject("pF[", f, "] is out of bounds at i =", i, ", id[i] =", id[i], ", trial[i] =", trial[i]);
        }
      }


      temp[4] = pF[1]; //Update pF[1] 
      temp[5] = pF[2]; //Update pF[2] 
      
      Prob_F[i, 1:5] = temp'; //Save current probabilities of sampling food f in Prob_F. ' transforms to row vector.
      
      //Up until sampling before this, after sampling "update" below this. 
      F[id[i], food_sample[i]] += 1; //increment the count of the sampled food type for the individual
      a_T_update = to_vector(F[id[i], 1:2]) + alpha_ind; //Update alpha by summing food type counts and their probability of being sampled (alpha)
      a_T[id[i]] = softmax(a_T_update); //Use softmax to convert to a valid probability vector
      //Check how this behaves. Does a_T reach estimated alpha for that session?
        
      //Check for NaN or inf values
      if (is_nan(lambda_ind) || is_inf(lambda_ind)) {
        reject("lambda_ind is invalid at i =", i, ", id[i] =", id[i], ", sess[i] =", sess[i]);
      }
      if (is_nan(pF[1]) || is_inf(pF[1]) || is_nan(pF[2]) || is_inf(pF[2])) {
        reject("pF is invalid at i =", i, ", id[i] =", id[i], ", sess[i] =", sess[i]);
      }
        
      //Here we start calculating stay prob
      //if(test[i] == 1){ //If this is a test trial
      //  A_row = to_vector(A[id[i], 1:2]); //Initialize vector of attraction scores of each id to each food type.    
      //  
      //  if (is_nan(A_row[1]) || is_inf(A_row[1]) || is_nan(A_row[2]) || is_inf(A_row[2])) {
      //    reject("A_row is invalid at i =", i, ", id[i] =", id[i], ", sess[i] =", sess[i]);
      //  }
        
      
      //  pF_x_A_row = pF .* A_row; //multiply the prob of finding a food type (pF) with the attraction score.
      //  scaled_pF_x_A_row = lambda_ind * pF_x_A_row; //Scale using individual lambda (I understand this as a value representing how closely an IDs choices align with underlying pronbabilities?)      
      //  pS = softmax(scaled_pF_x_A_row); //Transform to valid probability vector

      //  S_Prob[cum_trial[i]] = pS; //Save pS as S_Prob for each test trial

      //  temp[4] = pS[1];
      //  temp[5] = pS[2];

      //  Prob_S[cum_trial[i], 1:5] = temp';//Save Prob_S for each ID/test trial
        
      //}
      

      
      //Attraction updating (change this to update attraction based on whether or not food consumed.)
      //How does reinforcment learning affect switch decisions?
      //for(f in 1:N_food){ //For each food type f
      //  if(ate_sample[i] == 1){ //if the sampled food was eaten
      //    A[id[i], f] = (1 - phi_ind) * A[id[i], f] + phi_ind * Pay[id[i], f];//Increase Attraction for the sampled item 
      //    } else {
      //    A[id[i], f] = (1 - phi_ind) * A[id[i], f]; //Decrease Attraction for the sampled item (check what this line does exactly)
      //  }
      //}
      
    } 
  } 
} 
  
model{ 

  weight_p ~ exponential(1); //Weight parameter for food consumption             
  alpha_p ~ normal(0, 1); //Global intercept for food consumption                
  to_vector(V_p) ~ normal(0, 1); //Value parameters for each bird / food type      
  
  to_vector(logit_phi) ~  normal(0, 1); //attraction learning rate
  to_vector(log_lambda) ~  normal(0, 1); //scaling for prob
  
  for(j in 1:N_id){//for each ID
    for(s in 1:N_sess){//for each session
      raw_beta[j, s] ~ dirichlet(rep_vector(1, 2));//Flat prior for beta
      raw_alpha[j, s] ~ dirichlet(raw_beta[j, s]);//Alpha for each sess is tied to beta
    }
  }
  
  to_vector(z_ID) ~ normal(0, 1); 
  sigma_ID ~ exponential(1);      
  Rho_ID ~ lkj_corr_cholesky(4);  
//Likelihood functions
  //Predict food consumption for each ID and each food type
  for(j in 1:N_id){
    for(f in 1:N_food){
      total_eat[j, f] ~ poisson_log(alpha_p + log_dur[j] + weight_p * inv_logit(V_p[j, f]));
    }
  }
  
  
  //Stay / switch choice across all crows (for each test trial)
  //for(t in 1:T_test){
  //  choice[t] ~ categorical(S_Prob[t]);
  //}
  
}

generated quantities {

  matrix[N_id, N_sess] phi;     
  matrix[N_id, N_sess] lambda;  
  simplex[2] beta[N_id, N_sess]; 
  simplex[2] alpha[N_id, N_sess]; 

  for(j in 1:N_id){
    for(s in 1:N_sess){
      phi[j, s] = inv_logit(logit_phi[j, s]);
      lambda[j, s] = exp(log_lambda[j, s]);
      beta[j, s] = raw_beta[j, s];
      alpha[j, s] = raw_alpha[j, s];
    }
  }
}
