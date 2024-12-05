gen_OH_data <- function(N_id = 3, N_trials_sess = 10, N_sess = 10, N_trials_test_sess = 30, N_food = 2, group = "mixed"){
  #N_id <- 3 #Number of IDs
  #N_trials_sess <- 10 #Number of trials per session
  #N_sess <- 10 #Number of sessions (excluding test)
  #N_trials_test_sess <- 30 #Number of trials in the test session
  N_trials_id <- N_sess*N_trials_sess+N_trials_test_sess #Number of trials per ID (including test)
  N <- N_id*N_trials_id #Total number of trials
  #N_food <- 2 #Number of food types (here 2; low-value and high-value)
  T_test <- N_trials_test_sess*N_id #Total number of test trials 
  
  #total_eat is the total number of each food type consumed in the preference testing. high-value = 1
  total_eat <- matrix(0, nrow = N_id, ncol = N_food)
  for (i in 1:N_id) {  # Iterate over the IDs
    t <- table(factor(sample(c(1, 2), size = 100, replace = TRUE, prob = c(0.9, 0.1)), levels = 1:2)) #Skewed towards preference for high-value.
    total_eat[i, ] <- as.numeric(t)
  }
  
  #log_dur is the log10 total time spent in preference training for each crow in seconds (sum of seconds in each trial).
  log_dur <- matrix(0, nrow = N_id, ncol = 1)
  for (i in 1:N_id) {
    d <- sum(round(rnorm(100, mean = 10, sd = 5), 2))
    log_dur[i, ] <- log10(as.numeric(d))
  }
  
  id <- c()
  for(i in 1:N_id) id <- c(id, rep(i, N_trials_id))
  
  if (N_trials_test_sess == 0){#if there is no test session
    sess <- c()
    for(s in 1:N_sess) sess <- c(sess, rep(s, N_trials_sess))
    sess <- rep(sess, N_id)
    trial <- rep(rep(1:N_trials_sess, N_sess),N_id)
    cum_trial <- rep(1:(N_trials_id), N_id)

    if (group == "mixed"){
      n_switch <- 10
      
      food_sample <- c()
      for (i in 1:N_id){
        food_sample <- c(food_sample, sample(c(1,2), N_trials_id ,replace =TRUE, prob = c(0.5, 0.5)))
      }
    } else if (group == "uniform") {
      n_switch <- 3
      
      food_sample <- c()
      for (i in 1:N_id){
        food_sample <- c(food_sample, rep(c(rep(1,10), rep(2,10)),N_sess/2))
      }
    } else if (group == "both") {
      n_switch = 10
      food_sample_mixed <- c()
      for (i in 1:(round(N_id/2))){
        food_sample_mixed <- c(food_sample_mixed, sample(c(1,2), N_trials_id ,replace =TRUE, prob = c(0.5, 0.5)))
      }  
      n_switch = 3
      food_sample_uni <- c()
      for (i in 1:(N_id-round(N_id/2))){
        food_sample_uni <- c(food_sample_uni, rep(c(rep(1,10), rep(2,10)), N_sess/2))
      }     
      
      food_sample <- c(food_sample_mixed, food_sample_uni)
      
    } else {
      print("Group has to be specified as either mixed, uniform or both (for an equal mix of mixed and uniform).")
    } 
    
    test <- c(0)   
    choice <- c(0)
    
  } else { #if there is a test session
    sess <- c()
    for(s in 1:N_sess) sess <- c(sess, rep(s, N_trials_sess))
    sess <- rep(c(sess, rep(N_sess+1,N_trials_test_sess)), N_id)
    trial <- rep(c(rep(1:N_trials_sess, N_sess), 1:N_trials_test_sess),N_id)
    cum_trial <- rep(c(1:(N_trials_id-N_trials_test_sess), 1:N_trials_test_sess), N_id)
    
    if (group == "mixed"){
      n_switch <- 10
      
      food_sample <- c()
      choice <- c()
      for (i in 1:N_id){
        choice <- c(choice, c(sample(c(1,2), n_switch, replace = TRUE, prob = c(0.95, 0.05)), sample(c(1,2), N_trials_test_sess-n_switch, replace = TRUE, prob = c(0.5, 0.5))))
        food_sample <- c(food_sample,c(sample(c(1,2), N_trials_id-N_trials_test_sess ,replace =TRUE, prob = c(0.5, 0.5)), rep(2, N_trials_test_sess)))
      }
    } else if (group == "uniform") {
      n_switch <- 3
      
      food_sample <- c()
      choice <- c()
      for (i in 1:N_id){
        choice <- c(choice, c(sample(c(1,2), n_switch, replace = TRUE, prob = c(0.95, 0.05)), sample(c(1,2), N_trials_test_sess-n_switch, replace = TRUE, prob = c(0.5, 0.5))))
        food_sample <- c(food_sample,c(rep(c(rep(1,10), rep(2,10)),N_sess/2), rep(2, N_trials_test_sess)))
      }
    } else if (group == "both") {
      n_switch = 10
      choice_mixed <- c()
      food_sample_mixed <- c()
      for (i in 1:(round(N_id/2))){
        choice_mixed <- c(choice_mixed, c(sample(c(1,2), n_switch, replace = TRUE, prob = c(0.95, 0.05)), sample(c(1,2), N_trials_test_sess-n_switch, replace = TRUE, prob = c(0.5, 0.5))))
        food_sample_mixed <- c(food_sample_mixed,c(sample(c(1,2), N_trials_id-N_trials_test_sess ,replace =TRUE, prob = c(0.5, 0.5)), rep(2, N_trials_test_sess)))
      }  
      n_switch = 3
      choice_uniform <- c()
      food_sample_uni <- c()
      for (i in 1:(N_id-round(N_id/2))){
        choice_uniform <- c(choice_uniform, c(sample(c(1,2), n_switch, replace = TRUE, prob = c(0.95, 0.05)), sample(c(1,2), N_trials_test_sess-n_switch, replace = TRUE, prob = c(0.5, 0.5))))
        food_sample_uni <- c(food_sample_uni,c(rep(c(rep(1,10), rep(2,10)),N_sess/2), rep(2, N_trials_test_sess)))
      }     
      
      choice <- c(choice_mixed,choice_uniform)
      food_sample <- c(food_sample_mixed, food_sample_uni)
      
    } else {
      print("Group has to be specified as either mixed, uniform or both (for an equal mix of mixed and uniform).")
    } 
    
    test <- rep(c(rep(0, N_trials_id-N_trials_test_sess), rep(1, N_trials_test_sess)), N_id)
    N_sess <- N_sess+1 #Update N_sess to include the test (necessary for indexing in the model)
  }


  ate_sample <- c()
  for(i in food_sample){
    if(i == 1){ #If the food is high-value
      s <- sample(c(0,1), 1, prob = c(0.05, 0.95)) #95% chance to eat it
    } else if (i == 2) { #If the food is low value 
      s <- sample(c(0,1), 1, prob = c(0.9, 0.1)) #10% chance to eat it
    }
    ate_sample <- c(ate_sample, s)
  }
  
  
  dataList <- list(
    N = as.integer(N),
    N_id = as.integer(N_id),
    N_food = as.integer(N_food),
    N_sess = as.integer(N_sess),
    T_test = as.integer(T_test),
    
    total_eat = total_eat,
    log_dur = as.vector(log_dur),
    
    id = as.integer(id),
    sess = as.integer(sess),
    trial = as.integer(trial),
    cum_trial = as.integer(cum_trial),
    food_sample = as.integer(food_sample),
    ate_sample = as.integer(ate_sample),
    test = as.integer(test),
    choice = as.integer(choice)
  )
  
  return(dataList)
  
}
