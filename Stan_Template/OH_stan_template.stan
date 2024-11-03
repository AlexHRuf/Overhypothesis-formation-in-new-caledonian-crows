data{ 

  int N;           //Total number of observations       
  int N_id;        //Total number of birds      
  int N_food;      //Total number of food items (2)          
  int N_sess;      //Total number of sessions       
  int T_test;      //Total number of test trials across individuals        
  
  int total_eat[N_id, N_food]; //Total pieces of food item 1 and 2 eaten
  real log_dur[N_id];          //Total duration of preference test in seconds (log transformed)
  
  int id[N];        //Bird-specific individual identificaiton number (1 .. N)        
  int sess[N];      //Session identification number (1 .. N_sess)           
  int trial[N];     //Within session running trial account      
  int cum_trial[N]; //Running trial count that resets between phases      
  int sample[N];    //Food item sampled (1 = high, 2 = low)      
  int test[N];      //Whether the test session (0 = no, 1 = yes)   
  int choice[N];    //Whether stay or switch (1 = stay, 2 = switch) - note: cannot 0/1 dummy code dependent variabels in stan    

}
  
parameters{ 
    
  matrix[N_id, N_food] V_p;        
  real<lower=0> weight_p;          
  real alpha_p;                    
  
  matrix[N_id, N_sess] logit_phi;     
  matrix[N_id, N_sess] log_lambda;    
  simplex[2] raw_beta[N_id, N_sess];  
  simplex[2] raw_alpha[N_id, N_sess]; 
  
  matrix[3, N_id] z_ID;           
  vector<lower = 0>[3] sigma_ID;  
  cholesky_factor_corr[3] Rho_ID; 
 
} 

transformed parameters{ 

  matrix[N_id, 3] v_ID;      
  simplex[2] S_Prob[T_test]; 
  matrix[N, 5] alpha_T; 
  matrix[N, 5] Prob_F;  
  matrix[N, 5] Prob_S; 

  v_ID = (diag_pre_multiply(sigma_ID, Rho_ID) * z_ID)'; 

  { 

    matrix[N_id, N_food] F; 
    vector[2] a_T[N_id]; 
    matrix[N_id, N_food] Pay; 
    matrix[N_id, N_food] A;   
   
    for(j in 1:N_id){
    
      A[j, 1:2] = rep_vector(0.1, N_food)'; 
      
      for(f in 1:N_food){
        Pay[j, f] = inv_logit(V_p[j, f]); 
      }
    
    }

    for(i in 1:N){ 
    
      real phi_ind = inv_logit(logit_phi[id[i], sess[i]] + v_ID[id[i], 1]); 
      real lambda_ind = exp(log_lambda[id[i], sess[i]] + v_ID[id[i], 2]);  
     
      vector[2] alpha = raw_alpha[id[i], sess[i]]; 
      real alpha_unconstrained = log((alpha[1] + 1e-8) / (alpha[2] + 1e-8)); 
      real alpha_v_ID_sum = alpha_unconstrained + v_ID[id[i], 3]; 
      vector[2] alpha_v_ID_sum_transformed = softmax([alpha_v_ID_sum, 0]'); 
      vector[2] alpha_ind = alpha_v_ID_sum_transformed; 
      
      vector[5] temp;
      vector[2] pF; 
      vector[2] A_row;
      vector[2] pF_x_A_row;
      vector[2] scaled_pF_x_A_row;
      vector[2] pS;
      vector[2] a_T_update;
      
      if(trial[i] == 1){
        a_T[id[i]] = alpha_ind; 
        F[id[i], 1:2] = rep_vector(0, N_food)'; 
      }
      
      temp[1] = id[i];
      temp[2] = sess[i];
      temp[3] = trial[i];
      temp[4] = a_T[id[i]][1];
      temp[5] = a_T[id[i]][2];
      
      alpha_T[i, 1:5] = temp';
      
      for(f in 1:N_food){
        pF[f] = (F[id[i], f] + a_T[id[i]][f]) / (trial[i] + N_food);  
      }
      
      pF = pF / sum(pF); 

      temp[4] = pF[1];
      temp[5] = pF[2];
      
      Prob_F[i, 1:5] = temp';
      
      F[id[i], sample[i]] += 1; 
      a_T_update = to_vector(F[id[i], 1:2]) + alpha_ind; 
      a_T[id[i]] = softmax(a_T_update); 

      if(test[i] == 1){
        A_row = to_vector(A[id[i], 1:2]);        
        pF_x_A_row = pF .* A_row;                
        scaled_pF_x_A_row = lambda_ind * pF_x_A_row;        
        pS = softmax(scaled_pF_x_A_row);       

        S_Prob[cum_trial[i]] = pS; 

        temp[4] = pS[1];
        temp[5] = pS[2];

        Prob_S[i, 1:5] = temp';
        
      }
      
      for(f in 1:N_food){
        if(sample[i] == f){
          A[id[i], f] = (1 - phi_ind) * A[id[i], f] + phi_ind * Pay[id[i], f];
          } else {
          A[id[i], f] = (1 - phi_ind) * A[id[i], f];
        }
      }
      
    } 
  } 
} 
  
model{ 

  weight_p ~ exponential(1);             
  alpha_p ~ normal(0, 1);                
  to_vector(V_p) ~ normal(0, 1);         
  
  to_vector(logit_phi) ~  normal(0, 1);
  to_vector(log_lambda) ~  normal(0, 1);
  
  for(j in 1:N_id){
    for(s in 1:N_sess){
      raw_beta[j, s] ~ dirichlet(rep_vector(1, 2)); 
      raw_alpha[j, s] ~ dirichlet(raw_beta[j, s]);  
    }
  }
  
  to_vector(z_ID) ~ normal(0, 1); 
  sigma_ID ~ exponential(1);      
  Rho_ID ~ lkj_corr_cholesky(4);  

  for(j in 1:N_id){
    for(f in 1:N_food){
      total_eat[j] ~ poisson_log(alpha_p + log_dur[j] + weight_p * inv_logit(V_p[j, f]));
    }
  }

  for(t in 1:T_test){
    choice[t] ~ categorical(S_Prob[t]);
  }
  
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
