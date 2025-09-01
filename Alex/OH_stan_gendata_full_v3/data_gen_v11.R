# generate_stan_data.R  — matches the current Stan model
# ------------------------------------------------------------
# Simulates data for the model with:
#  • Scalar alpha per ID (alpha_s)
#  • beta_s simplex per ID, used as the starting β in test
#  • Monotone-to-endpoint updates for α and β (equal-weight schedules)
#  • Utilities: util = (exp(lambda_raw + u_lambda)) * U_board  -  exp(stickiness_raw) on switch
#  • Optional per-ID random effects for (lambda_raw, alpha_last_raw, beta_last_raw, phi_raw)
#  • Dynamic attraction update for LV as in the Stan model
#  • NEW: fixed attraction toggle (fixed_A) and attraction tracker (AF_t)
# ------------------------------------------------------------

generate_stan_data <- function(
    # --- sizes & indexing ---
  N_id = 16,
  N_board = 10,
  N_food = 2,              # must be 2
  seed = 12345,
  
  # --- sampling environment (matches your experimental design) ---
  total_samples = 120,     # total sampling pecks per ID across boards
  min_board     = 8,       # min pecks per board
  max_board     = 14,      # max pecks per board
  evidence_noise = c("multinomial","deterministic")[1],
  
  # --- initial OH for test (alpha_s scalar, beta_s simplex) ---
  mixed_alpha   = 8,       # α for "mixed" group
  uniform_alpha = 0.5,     # α for "uniform" group
  beta_init_mode = c("from_sampling","fixed"),  # how to set β_s
  fixed_beta      = c(0.5, 0.5),                # used if beta_init_mode="fixed"
  
  # --- test-phase parameter means (raw scales) ---
  lambda_raw     = 0.0,    # utility sensitivity (log-scale)
  phi_raw        = 0.0,    # attraction update rate (logit scale)
  stickiness_raw = 0.0,    # penalty (log-scale, exponentiated in utility)
  alpha_last_raw = 0.0,    # α endpoint raw (enter squared & negated)
  beta_last_raw  = 0.0,    # β HV-endpoint raw (enter squared & negated)
  
  # --- per-ID random effects for (λ, α_end, β_end, φ) ---
  heterogeneity = TRUE,
  sd_lambda = 0.3,
  sd_alpha_end = 0.3,
  sd_beta_end  = 0.3,
  sd_phi       = 0.3,
  corr_type = c("independent","exchangeable","manual"),
  rho_exchangeable = 0.0,                # one-number off-diagonal if exchangeable
  corr_matrix_manual = diag(4),          # 4x4 if corr_type="manual"
  
  # --- test-phase execution ---
  max_test_trials = 50,
  break_on_switch = F,
  
  # --- attractions (initial level & eating probability for test logging) ---
  A0 = c(0.9, 0.1),        # initial attractions (HV, LV), on 0..1
  fixed_A = NULL,          # NEW: if non-NULL (length-2), attractions are held fixed to this value
  p_eat_low  = 0.05,       # Pr(eat LV when seen on B11)
  p_eat_high = 0.95        # used for exposure tallies from sampling/PT
) {
  stopifnot(N_food == 2)
  set.seed(seed)
  beta_init_mode <- match.arg(beta_init_mode)
  corr_type      <- match.arg(corr_type)
  
  normalize <- function(v) { v <- as.numeric(v); s <- sum(v); if (!is.finite(s) || s <= 0) rep(1/length(v), length(v)) else v/s }
  softmax   <- function(x) { ex <- exp(x - max(x)); ex/sum(ex) }
  
  # ---- tiny rmvnorm via Cholesky ----
  rmvnorm_chol <- function(n, mu, Sigma) {
    L <- chol(Sigma)
    Z <- matrix(rnorm(n * length(mu)), nrow = n)
    sweep(Z %*% L, 2, mu, FUN = "+")
  }
  
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
  treatment <- rep(c("mixed","uniform"), length.out = N_id)
  
  env_prob <- vector("list", N_id)
  for (i in seq_len(N_id)) {
    if (treatment[i] == "mixed") {
      env_prob[[i]] <- replicate(N_board, c(0.5, 0.5), simplify = FALSE)
    } else {
      seq8   <- alternating_blocks()
      order10 <- c("H","L", seq8)  # H/L boards alternating blocks (no >2 in a row)
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
  
  # ---------- initial OH for test ----------
  # alpha_s: scalar per ID, from treatment
  alpha_s_vec <- ifelse(treatment == "mixed", mixed_alpha, uniform_alpha)
  # beta_s: simplex per ID
  if (beta_init_mode == "from_sampling") {
    beta_s_mat <- t(apply(y_sampling, 1, function(a) normalize(colSums(a))))
  } else {
    stopifnot(length(fixed_beta) == 2, all(fixed_beta >= 0))
    beta_s_mat <- matrix(normalize(fixed_beta), nrow = N_id, ncol = 2, byrow = TRUE)
  }
  
  # ---------- per-ID random effects for (λ, α_end, β_end, φ) ----------
  mu4   <- c(lambda_raw, alpha_last_raw, beta_last_raw, phi_raw)
  sds4  <- c(sd_lambda,  sd_alpha_end,  sd_beta_end,  sd_phi)
  if (corr_type == "independent") {
    Corr4 <- diag(4)
  } else if (corr_type == "exchangeable") {
    rho <- rho_exchangeable
    Corr4 <- matrix(rho, 4, 4); diag(Corr4) <- 1
  } else {
    Corr4 <- corr_matrix_manual
  }
  Sigma4 <- diag(sds4) %*% Corr4 %*% diag(sds4)
  
  # *** FIXED: robust construction when heterogeneity == FALSE ***
  if (heterogeneity) {
    RE <- rmvnorm_chol(N_id, mu4, Sigma4)
  } else {
    RE <- matrix(mu4, nrow = N_id, ncol = 4, byrow = TRUE)  # each row = mu4
  }
  colnames(RE) <- c("lambda_raw_i","alpha_last_raw_i","beta_last_raw_i","phi_raw_i")
  
  # Natural scales used in utilities and attraction update
  if (heterogeneity) {
    L_id <- exp(RE[, "lambda_raw_i"])     # length N_id
    P_id <- plogis(RE[, "phi_raw_i"])     # length N_id
  } else {
    L_id <- rep(exp(lambda_raw), N_id)    # repeat same value for all IDs
    P_id <- rep(plogis(phi_raw), N_id)
  }
  S <- exp(stickiness_raw) 
  
  # ---------- trackers ----------
  P_S_tracker       <- matrix(NA_real_, nrow = N_id, ncol = max_test_trials)
  theta_b11_tracker <- array(NA_real_, dim = c(N_id, max_test_trials, N_food))
  theta_b12_tracker <- array(NA_real_, dim = c(N_id, max_test_trials, N_food))
  alpha_tracker     <- matrix(NA_real_, nrow = N_id, ncol = max_test_trials)    # scalar α_t
  beta_tracker      <- array(NA_real_, dim = c(N_id, max_test_trials, N_food))  # β_t
  AF_tracker        <- array(NA_real_, dim = c(N_id, max_test_trials, N_food))  # NEW: attractions
  
  # ---------- outputs (Stan data fields) ----------
  Switch      <- integer(0)  # 1=stay, 2=switch
  id_test     <- integer(0)
  trial_num   <- integer(0)
  y_test      <- integer(0)  # food observed on B11 (LV=2) until switch
  y_test_ate  <- integer(0)  # ate? 0/1
  T_id        <- integer(N_id)
  
  # ---------- test phase ----------
  for (i in seq_len(N_id)) {
    # per-ID initial state
    a_first <- alpha_s_vec[i]               # scalar α
    b_first <- as.numeric(beta_s_mat[i, ])  # length-2 simplex (HV, LV)
    
    # *** FIXED: endpoints use RE only if heterogeneity is TRUE ***
    alpha_end_raw_i <- if (heterogeneity) RE[i, "alpha_last_raw_i"] else alpha_last_raw
    beta_end_raw_i  <- if (heterogeneity) RE[i,  "beta_last_raw_i"] else beta_last_raw
    
    # per-ID endpoints using squared-shrink parameterization
    a_last     <- a_first    * exp(- alpha_end_raw_i^2)
    p_HV_last  <- b_first[1] * exp(- beta_end_raw_i^2)
    b_last     <- c(p_HV_last, 1 - p_HV_last)
    
    # equal-weight schedules across planned horizon (simple monotone schedules)
    T_plan <- max_test_trials
    wA <- rep(1 / T_plan, T_plan)
    wB <- rep(1 / T_plan, T_plan)
    
    counts_b11 <- c(0.0, 0.0)
    # NEW: fixed attraction toggle
    if (is.null(fixed_A)) {
      A_F <- as.numeric(A0)                 # start from A0, allow learning
      do_learn_A <- TRUE
    } else {
      stopifnot(length(fixed_A) == 2, all(is.finite(fixed_A)), all(fixed_A >= 0), all(fixed_A <= 1))
      A_F <- as.numeric(fixed_A)            # hold attractions fixed
      do_learn_A <- FALSE
    }
    t_done <- 0
    
    for (t in 1:max_test_trials) {
      # cum weights use *previous* trials so trial 1 = start
      cumA <- if (t >= 2) sum(wA[1:(t-1)]) else 0
      cumB <- if (t >= 2) sum(wB[1:(t-1)]) else 0
      
      # current α_t, β_t via interpolation to endpoints
      alpha_t <- a_first + (a_last - a_first) * cumA
      beta_t  <- b_first + (b_last - b_first) * cumB
      
      # pseudo-counts and predictive per board
      OH <- as.numeric(alpha_t) * beta_t        # vector length 2
      theta_b11 <- OH + counts_b11
      theta_b12 <- OH
      P_F_b11 <- normalize(theta_b11)
      P_F_b12 <- normalize(theta_b12)
      
      # utilities → stay/switch probabilities
      U_b11 <- sum(P_F_b11 * A_F)
      U_b12 <- sum(P_F_b12 * A_F)
      util  <- c(L_id[i] * U_b11,
                 L_id[i] * U_b12 - S)
      P     <- softmax(util)  # c(P_stay, P_switch)
      
      # record a test observation BEFORE outcome
      y_test     <- c(y_test, 2L)                   # LV on board 11
      ate        <- rbinom(1, 1, p_eat_low)
      y_test_ate <- c(y_test_ate, ate)
      
      # trackers (pre-update, consistent with existing trackers)
      theta_b11_tracker[i, t, ] <- P_F_b11
      theta_b12_tracker[i, t, ] <- P_F_b12
      alpha_tracker[i, t]       <- alpha_t
      beta_tracker[i, t, ]      <- beta_t
      AF_tracker[i, t, ]        <- A_F        # NEW: attractions snapshot
      P_S_tracker[i, t]         <- P[2]
      
      # realised choice
      choice    <- sample(1:2, 1, prob = P)
      Switch    <- c(Switch, choice)
      id_test   <- c(id_test, i)
      trial_num <- c(trial_num, t)
      
      # stop after switch (like Stan data collection: no more B11 trials)
      if (break_on_switch && choice == 2) {
        t_done <- t
        break
      }
      
      # stayed → observe LV on b11, update counts
      counts_b11[2] <- counts_b11[2] + 1.0
      
      # update attractions for LV exactly like the model (if not fixed)
      if (do_learn_A) {
        P_upd <- P_id[i]   # inv_logit(phi_raw + u_phi)
        if (ate == 1) {
          A_F[2] <- (1 - P_upd) * A_F[2] + P_upd * 0.25
        } else {
          A_F[2] <- (1 - P_upd) * A_F[2]  # no reward
        }
        # (HV unchanged in the Stan model)
      }
    }
    
    if (break_on_switch) {
      if (t_done == 0) t_done <- max_test_trials
      T_id[i] <- t_done
    } else {
      T_id[i] <- max_test_trials
    }
  }
  
  T_max <- max(T_id)
  
  # ---------- Preference test & totals (exposure-style) ----------
  PT_log_dur <- log(runif(N_id, 300, 600))
  PT_trials  <- sample(60:80, N_id, TRUE)
  PT_ate_HV  <- rbinom(n = N_id, size = PT_trials, prob = p_eat_high)
  PT_ate_LV  <- rbinom(n = N_id, size = PT_trials, prob = p_eat_low)
  PT_total_eat <- cbind(PT_ate_HV, PT_ate_LV)
  
  # Totals from sampling (exposure), shape: N_id x 2 (HV, LV)
  samp_counts <- apply(y_sampling, c(1,3), sum)
  samp_ate_HV <- rbinom(n = N_id, size = samp_counts[, 1], prob = p_eat_high)
  samp_ate_LV <- rbinom(n = N_id, size = samp_counts[, 2], prob = p_eat_low)
  samp_ate <- cbind(samp_ate_HV, samp_ate_LV)
  
  ALL_total_eat <- PT_total_eat + samp_ate
  ALL_log_dur   <- PT_log_dur * 2
  
  # ensure integer storage for Stan
  storage.mode(PT_total_eat) <- "integer"
  storage.mode(ALL_total_eat) <- "integer"
  storage.mode(y_sampling)    <- "integer"
  
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
    
    # extras (useful for checks)
    treatment = treatment,
    alpha_s   = alpha_s_vec,
    beta_s    = beta_s_mat,
    RE_draws  = RE,                # per-ID draws for (λ, α_end, β_end, φ)
    L_id      = L_id,
    stickiness_raw = stickiness_raw,
    P_update  = P_id,
    trackers = list(
      P_S = P_S_tracker,
      theta_b11 = theta_b11_tracker,
      theta_b12 = theta_b12_tracker,
      alpha_t = alpha_tracker,
      beta_t  = beta_tracker,
      AF_t    = AF_tracker        # NEW: attractions tracker
    ),
    # params used
    args = list(
      total_samples = total_samples,
      min_board = min_board,
      max_board = max_board,
      evidence_noise = evidence_noise,
      mixed_alpha = mixed_alpha,
      uniform_alpha = uniform_alpha,
      beta_init_mode = beta_init_mode,
      fixed_beta = fixed_beta,
      lambda_raw = lambda_raw,
      phi_raw = phi_raw,
      stickiness_raw = stickiness_raw,
      alpha_last_raw = alpha_last_raw,
      beta_last_raw  = beta_last_raw,
      heterogeneity = heterogeneity,
      sd_lambda = sd_lambda,
      sd_alpha_end = sd_alpha_end,
      sd_beta_end  = sd_beta_end,
      sd_phi       = sd_phi,
      corr_type = corr_type,
      rho_exchangeable = rho_exchangeable,
      corr_matrix_manual = corr_matrix_manual,
      max_test_trials = max_test_trials,
      break_on_switch = break_on_switch,
      A0 = A0,
      fixed_A = fixed_A,           # NEW: echo choice
      p_eat_low = p_eat_low,
      p_eat_high = p_eat_high
    )
  )
}

# --- Example call (kept as in your message) ---
dataList <- generate_stan_data(
  # sizes & indexing
  N_id = 30,
  N_board = 10,
  N_food = 2,
  seed = 264346,
  
  # sampling environment
  total_samples = 120,
  min_board = 8,
  max_board = 14,
  evidence_noise = "multinomial",   # or "deterministic"
  
  # initial OH for test
  mixed_alpha   = 10,                # α for "mixed" group
  uniform_alpha = 0.5,              # α for "uniform" group
  beta_init_mode = "fixed", # or "from_sampling"
  fixed_beta      = c(0.5, 0.5),    # only used if beta_init_mode = "fixed"
  
  # test-phase parameter means (raw scales)
  lambda_raw     = 3,   # log-sensitivity (L = exp(lambda_raw + u_lambda))
  phi_raw        = 1,   # attraction update (P = inv_logit(phi_raw + u_phi))
  stickiness_raw = 1.5, # switch penalty (S = exp(stickiness_raw))
  alpha_last_raw = 1,  # α endpoint (enter squared & negated in sim)
  beta_last_raw  = 1,  # β HV endpoint (enter squared & negated in sim)
  
  # per-ID random effects (λ, α_end, β_end, φ)
  heterogeneity = TRUE,
  sd_lambda   = 0.2,
  sd_alpha_end = 0.30,
  sd_beta_end  = 0.30,
  sd_phi       = 0.50,
  corr_type = "exchangeable",       # "independent", "exchangeable", or "manual"
  rho_exchangeable = 0.2,           # off-diagonal correlation if exchangeable
  corr_matrix_manual = diag(4),     # used only if corr_type = "manual"
  
  # test-phase execution
  max_test_trials = 50,
  break_on_switch = TRUE,
  
  # attractions & eating
  A0 = c(0.8, 0.2),        # initial attractions (HV, LV)
  fixed_A = NULL,
  p_eat_low  = 0.05,
  p_eat_high = 0.95
)

# Quick checks (unchanged)

df <- tibble(T_id = dataList$T_id, treatment = dataList$treatment)
ggplot(df, aes(y=T_id, fill = treatment))+geom_boxplot()+theme_minimal()

plot_tracker(dataList, dataList$trackers$alpha_t, y_lims = c(0,10))
plot_tracker(dataList, dataList$trackers$beta_t, food = 2)
plot_tracker(dataList, dataList$trackers$theta_b11, food = 2)
plot_tracker(dataList, dataList$trackers$theta_b12, food = 2)
plot_tracker(dataList, dataList$trackers$AF_t, food = 1)
plot_tracker(dataList, dataList$trackers$AF_t, food = 2)
plot_tracker(dataList, dataList$trackers$P_S)

