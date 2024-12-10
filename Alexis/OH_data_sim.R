generate_stan_data <- function(N_id = 16) {
  
  #Parameters
  N_food <- 2 # Number of food items (high = 1, low = 2)
  N_treat <- 2 # Number of treatments (mixed = 1, uniform = 2)
  N_sess <- 11 # Total number of sessions
  
  #Initialize vectors
  id <- integer(0)
  treat <- integer(0)
  sess <- integer(0)
  trial <- integer(0)
  sample_vec <- integer(0)
  choice <- integer(0)
  test <- integer(0)
  pseudo_row <- integer(0)
  phase <- integer(0)
  
  cum_trial_phase1 <- 0 #Initialize cumulative trial counter for phase 1
  cum_trial_phase2 <- 0 #Initialize cumulative trial counter for phase 2
  
  for(bird in 1:N_id){
    #Assign treatment
    bird_treat <- ifelse(bird <= N_id / 2, 1, 2) #Half in treatment 1, half in treatment 2
    
    #Allocate trials for sessions 1-10
    trials_per_session <- sample(1:20, N_sess - 1, replace = TRUE)
    trials_per_session <- round(trials_per_session / sum(trials_per_session) * 100)
    trials_per_session[length(trials_per_session)] <- 100 - sum(trials_per_session[-length(trials_per_session)])
    
    #Determine balanced session assignments for uniform treatment
    if (bird_treat == 2){
      session_assignments <- c(rep(1, (N_sess - 1) / 2), rep(2, (N_sess - 1) / 2)) #Equal high (1) and low (2)
      session_assignments <- sample(session_assignments) #Shuffle to randomize
      while (any(rle(session_assignments)$lengths > 2)){
        session_assignments <- sample(session_assignments) #Ensure no more than 2 consecutive repeats
      }
    }
    
    for(s in 1:(N_sess - 1)){
      id <- c(id, rep(bird, trials_per_session[s]))
      treat <- c(treat, rep(bird_treat, trials_per_session[s]))
      sess <- c(sess, rep(s, trials_per_session[s]))
      trial <- c(trial, 1:trials_per_session[s])
      test <- c(test, rep(0, trials_per_session[s])) #Test = 0 for sessions 1-10
      phase <- c(phase, rep(1, trials_per_session[s])) #Phase = 1 for sessions 1-10
      pseudo_row <- c(pseudo_row, cum_trial_phase1 + (1:trials_per_session[s]))
      cum_trial_phase1 <- cum_trial_phase1 + trials_per_session[s]
      
      if (bird_treat == 1){
        #Mixed treatment: Random sampling
        sample_vec <- c(sample_vec, sample(1:N_food, trials_per_session[s], replace = TRUE))
      } else {
        #Uniform treatment: Randomized balanced high/low assignments
        food_choice <- session_assignments[s]
        sample_vec <- c(sample_vec, rep(food_choice, trials_per_session[s]))
      }
    }
    
    #Add session 11 (test phase)
    session_11_trials <- if(bird_treat == 1){
      sample(15:20, 1)
    } else {
      sample(3:8, 1)
    }
    
    id <- c(id, rep(bird, session_11_trials + 1)) #+1 for the final switch trial
    treat <- c(treat, rep(bird_treat, session_11_trials + 1))
    sess <- c(sess, rep(N_sess, session_11_trials + 1))
    trial <- c(trial, 1:(session_11_trials + 1))
    test <- c(test, rep(1, session_11_trials + 1)) #Test = 1 for session 11
    phase <- c(phase, rep(2, session_11_trials + 1)) #Phase = 2 for session 11
    pseudo_row <- c(pseudo_row, cum_trial_phase2 + (1:(session_11_trials + 1)))
    cum_trial_phase2 <- cum_trial_phase2 + session_11_trials + 1
    
    sample_vec <- c(sample_vec, rep(2, session_11_trials), 1) #Session 11 all low quality, then switch
    choice <- c(choice, rep(1, session_11_trials), 2) #Stay = 1, Switch = 2
  }
  
  #Combine into a data frame
  data_frame <- data.frame(
    id = id,
    treat = treat,
    sess = sess,
    trial = trial,
    sample = sample_vec,
    test = test,
    phase = phase,
    pseudo_row = pseudo_row
  )
  
  #Create additional list elements
  total_eat <- matrix(c(sample(75:90, N_id, replace = TRUE), sample(0:10, N_id, replace = TRUE)), 
                      nrow = N_id, ncol = N_food, byrow = FALSE)
  log_dur <- log(runif(N_id, min = 300, max = 600))
  
  #Return list with data frame and additional elements
  return(list(
    data_frame = data_frame,
    total_eat = total_eat,
    log_dur = log_dur,
    choice = choice
  ))
}

#Run
data_sim <- generate_stan_data(N_id = 10)