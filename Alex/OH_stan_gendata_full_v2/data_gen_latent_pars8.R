generate_stan_data <- function(
    N_id = 16,
    N_board = 10,
    N_food = 2,
    seed = 12345,
    
    # --- OH initialisation for test (stand-ins for S_alpha, S_beta) ---
    mixed_alpha_total   = 30,
    uniform_alpha_total = 0.5,
    fixed_global_counts = c(60, 60),   # if NULL, per-ID totals from sampling are used
    
    # --- Sensitivity (lambda) and OH learning/smoothing ---
    lambda_raw = 0.0,   # L = exp(lambda_raw)
    rho_a = 0.20,       # used in dirichlet alpha mode
    rho_b = 0.20,       # smoothing when beta follows alpha
    
    # --- Fixed attractions on outcome scale (no learning here) ---
    fixed_A = c(0.9, 0.1),
    
    # --- Trial effect & stickiness (asymmetric; added to "switch" only) ---
    trial_effect = 0.0,   # (kept exactly as in your working version) *t*
    stickiness   = 0.0,
    
    # --- Update toggles ---
    alpha_mode = c("dirichlet","decay"),
    alpha_decay = 0.05,
    alpha_baseline = c(1, 1),
    
    beta_mode = c("alpha","global_counts"),
    
    # --- Sampling phase controls ---
    total_samples = 120,
    min_board     = 8,
    max_board     = 14,
    evidence_noise = c("multinomial", "deterministic")[1],
    
    # --- Test phase ---
    max_test_trials = 50,
    break_on_switch = TRUE,
    
    # --- Test-phase eating (for y_test_ate; no attraction learning) ---
    p_eat_low  = 0.05,
    p_eat_high = 0.95   # present for completeness; not used in this setup
){
  stopifnot(N_food == 2)
  set.seed(seed)
  
  alpha_mode <- match.arg(alpha_mode)
  beta_mode  <- match.arg(beta_mode)
  evidence_noise <- match.arg(evidence_noise)
  
  normalize <- function(v) { v <- as.numeric(v); s <- sum(v); if (!is.finite(s) || s <= 0) rep(1/length(v), length(v)) else v/s }
  softmax   <- function(x) { ex <- exp(x - max(x)); ex/sum(ex) }
  
  # ---------- helpers for sampling phase ----------
  random_session_counts <- function(N_board, min_board, max_board, total) {
    base <- rep(min_board, N_board)
    rem  <- total - sum(base)
    if (rem < 0) stop("total_samples too small for min_board*N_board")
    caps <- rep(max_board - min_board, N_board)
    while (rem > 0) {
      i <- sample.int(N_board, 1)
      if (caps[i] > 0) { base[i] <- base[i] + 1; caps[i] <- caps[i] - 1; rem <- rem - 1 }
    }
    base
  }
  alternating_blocks <- function() {
    v <- c(rep("H",4), rep("L",4))
    tries <- 0
    repeat {
      w <- sample(v, length(v), replace = FALSE)
      if (all(rle(w)$lengths <= 2)) return(w)
      tries <- tries + 1
      if (tries > 5000) stop("Failed to make uniform sequence with <=2 in a row.")
    }
  }
  
  # ---------- treatment & sampling env ----------
  treatment <- rep(c("mixed","uniform"), each = N_id/2)
  
  env_prob <- vector("list", N_id)
  for (i in seq_len(N_id)) {
    if (treatment[i] == "mixed") {
      env_prob[[i]] <- replicate(N_board, c(0.5, 0.5), simplify = FALSE)
    } else {
      seq8 <- alternating_blocks()
      order10 <- c("H","L", seq8)
      env_prob[[i]] <- lapply(order10, function(tag) if (tag=="H") c(1,0) else c(0,1))
    }
  }
  
  samples_per_board <- t(replicate(
    N_id, random_session_counts(N_board, min_board, max_board, total_samples)
  ))
  
  y_sampling <- array(0L, dim = c(N_id, N_board, N_food))
  for (i in seq_len(N_id)) for (b in seq_len(N_board)) {
    size <- samples_per_board[i,b]; p <- env_prob[[i]][[b]]
    if (evidence_noise == "multinomial") {
      k1 <- rbinom(1, size, p[1]); y_sampling[i,b,1] <- k1; y_sampling[i,b,2] <- size-k1
    } else {
      y_sampling[i,b,] <- as.integer(round(p * size))
      diff <- size - sum(y_sampling[i,b,])
      if (diff != 0) { j <- which.max(p); y_sampling[i,b,j] <- y_sampling[i,b,j] + diff }
    }
  }
  
  # ---------- initial OH state for test (per ID) ----------
  alpha0 <- matrix(NA_real_, nrow = N_id, ncol = N_food)
  beta0  <- matrix(NA_real_, nrow = N_id, ncol = N_food)
  for (i in seq_len(N_id)) {
    alpha0[i, ] <- if (treatment[i] == "mixed")
      c(mixed_alpha_total, mixed_alpha_total) else c(uniform_alpha_total, uniform_alpha_total)
    beta0[i, ] <- if (is.null(fixed_global_counts)) normalize(colSums(y_sampling[i, , ])) else fixed_global_counts / sum(fixed_global_counts)
  }
  
  # ---------- constants for test ----------
  L <- exp(lambda_raw)
  A_F <- as.numeric(fixed_A)   # fixed attractions (no learning)
  
  # ---------- trackers (diagnostics) ----------
  P_S_tracker       <- matrix(NA_real_, nrow = N_id, ncol = max_test_trials)      # P(switch)
  theta_b11_tracker <- array(NA_real_, dim = c(N_id, max_test_trials, N_food))     # P_F_b11
  theta_b12_tracker <- array(NA_real_, dim = c(N_id, max_test_trials, N_food))     # P_F_b12
  alpha_tracker     <- array(NA_real_, dim = c(N_id, max_test_trials, N_food))
  beta_tracker      <- array(NA_real_, dim = c(N_id, max_test_trials, N_food))
  
  # ---------- outputs (Stan data fields) ----------
  Switch      <- integer(0)  # 1=stay, 2=switch
  id_test     <- integer(0)
  trial_num   <- integer(0)
  y_test      <- integer(0)  # food observed on b11 (LV=2) until switch
  y_test_ate  <- integer(0)  # ate? 0/1
  T_id        <- integer(N_id)
  
  # ---------- test phase ----------
  for (i in seq_len(N_id)) {
    alpha_s <- as.numeric(alpha0[i, ])
    beta_s  <- as.numeric(beta0[i, ])
    counts_b11 <- c(0.0, 0.0)
    global_counts <- if (is.null(fixed_global_counts)) colSums(y_sampling[i, , ]) else as.numeric(fixed_global_counts)
    t_done <- 0
    
    for (t in 1:max_test_trials) {
      # thetas
      theta_b11  <- alpha_s * beta_s + counts_b11
      theta_b12  <- alpha_s * beta_s
      
      # P_F from theta (matches Stan)
      P_F_b11 <- normalize(theta_b11)
      P_F_b12 <- normalize(theta_b12)
      
      # utilities with trial effect ONLY on "switch" (kept exactly as your working version)
      U_b11 <- sum(P_F_b11 * A_F)
      U_b12 <- sum(P_F_b12 * A_F)
      util  <- c(L * U_b11, L * U_b12 + trial_effect * t - stickiness)
      P     <- softmax(util)  # c(P_stay, P_switch)
      
      # record a test observation BEFORE we know switch outcome:
      # b11 shows LV piece (2); we log y_test and whether it was eaten
      y_test     <- c(y_test, 2L)
      ate        <- rbinom(1, 1, p_eat_low)
      y_test_ate <- c(y_test_ate, ate)
      
      # track (pre-choice)
      theta_b11_tracker[i, t, ] <- P_F_b11
      theta_b12_tracker[i, t, ] <- P_F_b12
      alpha_tracker[i, t, ]     <- alpha_s
      beta_tracker[i, t, ]      <- beta_s
      P_S_tracker[i, t]         <- P[2]
      
      # realised choice
      choice <- sample(1:2, 1, prob = P)
      Switch    <- c(Switch, choice)
      id_test   <- c(id_test, i)
      trial_num <- c(trial_num, t)
      
      # if switch: stop after logging this trial's b11 LV (no further updates)
      if (break_on_switch && choice == 2) { 
        # note: on switch trial we recorded LV with ate above; that's fine (no learning from attraction)
        t_done <- t
        break 
      }
      
      # stay → observe LV on b11, update counts and OH
      counts_b11[2] <- counts_b11[2] + 1.0
      
      # --- ALPHA UPDATE ---
      if (alpha_mode == "decay") {
        alpha_s <- (1 - alpha_decay) * alpha_s + alpha_decay * as.numeric(alpha_baseline)
      } else { # "dirichlet"
        obs <- c(0.0, 1.0)                # saw LV
        alpha_new <- alpha_s + obs
        alpha_s   <- (1 - rho_a) * alpha_s + rho_a * alpha_new
      }
      
      # --- BETA UPDATE ---
      if (beta_mode == "alpha") {
        beta_new <- alpha_s / sum(alpha_s)
      } else { # "global_counts"
        global_counts <- global_counts + c(0, 1)  # increment LV count after stay
        beta_new <- normalize(global_counts)
      }
      beta_s <- (1 - rho_b) * beta_s + rho_b * beta_new
    }
    
    # Trials completed for this bird
    if (break_on_switch) {
      if (t_done == 0) t_done <- max_test_trials
      T_id[i] <- t_done
    } else {
      T_id[i] <- max_test_trials
    }
  }
  
  T_max <- max(T_id)
  
  # ---------- Preference test & totals (simple, fixed-attraction-compatible) ----------
  # We keep this simple to satisfy the Stan data{} interface, consistent with your earlier working approach.
  PT_log_dur <- log(runif(N_id, 300, 600))
  # Heavily prefer HV in PT:
  PT_total_eat <- cbind(sample(75:90, N_id, TRUE), sample(0:10, N_id, TRUE))
  storage.mode(PT_total_eat) <- "integer"
  
  # Evidence totals from sampling (per-ID, per-food)
  samp_counts <- apply(y_sampling, c(1,3), sum)
  # Add them to PT to form ALL_* (kept simple and consistent with previous code)
  ALL_total_eat <- PT_total_eat + samp_counts
  storage.mode(ALL_total_eat) <- "integer"
  ALL_log_dur <- PT_log_dur * 2
  
  # ---------- Final assembled outputs for Stan ----------
  list(
    # sizes / indexing
    N_test = length(y_test),
    N_food = N_food,
    N_id   = N_id,
    N_board = N_board,
    T_id = as.integer(T_id),
    T_max = as.integer(T_max),
    
    # preference test
    PT_total_eat = PT_total_eat,
    PT_log_dur   = PT_log_dur,
    
    # sampling phase
    y_sampling = y_sampling,
    
    # pref + sampling totals
    ALL_total_eat = ALL_total_eat,
    ALL_log_dur   = ALL_log_dur,
    
    # test phase
    y_test      = as.integer(y_test),
    y_test_ate  = as.integer(y_test_ate),
    id          = as.integer(id_test),
    trial       = as.integer(trial_num),
    Switch      = as.integer(Switch),
    
    # extras / trackers (handy but not used by Stan)
    treatment = treatment,
    alpha_init = alpha0,
    beta_init  = beta0,
    P_S_tracker = P_S_tracker,
    theta_b11_tracker = theta_b11_tracker,
    theta_b12_tracker = theta_b12_tracker,
    alpha_tracker = alpha_tracker,
    beta_tracker  = beta_tracker,
    
    # parameters used (recorded for reproducibility)
    lambda_raw = lambda_raw,
    rho_a = rho_a, rho_b = rho_b,
    fixed_A = fixed_A,
    trial_effect = trial_effect,
    stickiness = stickiness,
    alpha_mode = alpha_mode,
    alpha_decay = alpha_decay,
    alpha_baseline = alpha_baseline,
    beta_mode = beta_mode,
    break_on_switch = break_on_switch
  )
}

dataList <- generate_stan_data(N_id = 16 , 
                                 mixed_alpha_total = 20, 
                                 uniform_alpha_total = 1,
                                 fixed_global_counts = NULL,
                                 alpha_mode = "dirichlet", #"dirichlet" or "decay"
                                 alpha_decay = 0.1, #if alpha_mode = "decay"
                                 beta_mode = "alpha", #"global_counts" or "alpha"
                                 rho_a = 1, #if alpha_mode = "dirichlet"
                                 rho_b = 1, #if beta_mode = "alpha"
                                 lambda_raw = 3,
                                 trial_effect = -0.2,
                                 stickiness = 0,
                                 break_on_switch = F
)

source("plot_tracker.R")
dataList$T_id

plot_tracker(dataList = dataList,tracker = "P_S_tracker")
plot_tracker(dataList = dataList,food = 2,tracker = "theta_b11_tracker")
plot_tracker(dataList = dataList,food = 2,tracker = "theta_b12_tracker")

plot_tracker(dataList = dataList,food = 2,tracker = "alpha_tracker", y_lims = c(0,50))
plot_tracker(dataList = dataList,food = 2,tracker = "beta_tracker")

dataList$T_id

