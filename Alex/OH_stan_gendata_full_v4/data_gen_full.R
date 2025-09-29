#data_gen_full.R
generate_stan_data_full <- function(
    # --- sizes & indexing ---
  N_id = 16,
  N_board = 10,
  N_food = 2,              # must be 2
  seed = 12345,
  
  # --- sampling environment ---
  total_samples = 120,
  min_board     = 8,
  max_board     = 14,
  
  # --- initial OH for test (alpha_s scalar, beta_s simplex) ---
  alpha_mixed   = 8,       # α for "mixed" group (natural)
  alpha_uniform = 0.5,     # α for "uniform" group (natural)
  beta_mixed    = c(0.5, 0.5),  # β for "mixed" group (simplex)
  beta_uniform  = c(0.5, 0.5),  # β for "uniform" group (simplex)
  
  # --- treatment-level means (OUTCOME / NATURAL scales) ---
  # λ (utility sensitivity): outcome L > 0
  lambda_mixed   = 1.0,
  lambda_uniform = 1.0,
  # φ (attraction update):  outcome P in (0,1)
  phi_mixed      = 0.5,
  phi_uniform    = 0.5,
  # κ (stickiness / switch cost): outcome K >= 0
  kappa_mixed    = 1.0,
  kappa_uniform  = 1.0,
  
  # --- per-ID random effects sd's for (λ, φ, κ) ON RAW SCALE ---
  heterogeneity = TRUE,
  sd_lambda = 0.3,
  sd_phi    = 0.3,
  sd_kappa  = 0.3,
  corr_type = c("independent","exchangeable","manual"),
  rho_exchangeable = 0.0,
  corr_matrix_manual = diag(3),   # 3x3 if corr_type="manual"
  
  # --- test-phase execution ---
  max_test_trials = 50,
  break_on_switch = FALSE,
  
  # --- attractions & eating ---
  A0 = c(0.9, 0.1),        # initial attractions (HV, LV), on 0..1
  fixed_A = NULL,          # if non-NULL, attractions held fixed (φ has no effect)
  p_eat_low  = 0.05,
  p_eat_high = 0.95,
  
  # --- exponential decays (matches Stan) ---
  d_alpha = 0.01,          # 0..0.2 like Stan priors’ support
  d_beta  = 0.01
) {
  stopifnot(N_food == 2)
  set.seed(seed)
  corr_type <- match.arg(corr_type)
  
  normalize <- function(v) {
    v <- as.numeric(v); s <- sum(v)
    if (!is.finite(s) || s <= 0) rep(1/length(v), length(v)) else v/s
  }
  softmax   <- function(x) { ex <- exp(x - max(x)); ex/sum(ex) }
  invlogit  <- function(x) 1/(1+exp(-x))
  logit <- function(p) log(p/(1-p))
  
  # --- sanity checks on outcome-scale inputs ---
  if (lambda_mixed <= 0 || lambda_uniform <= 0)
    stop("lambda_* must be > 0 (outcome scale).")
  if (phi_mixed <= 0 || phi_mixed >= 1 || phi_uniform <= 0 || phi_uniform >= 1)
    stop("phi_* must be in (0,1) (outcome scale).")
  if (kappa_mixed < 0 || kappa_uniform < 0)
    stop("kappa_* must be >= 0 (outcome scale).")
  
  # ---- tiny rmvnorm via Cholesky with per-ID means ----
  rmvnorm_chol_means <- function(mu_mat, Sigma) {
    N <- nrow(mu_mat); K <- ncol(mu_mat)
    L <- chol(Sigma)
    Z <- matrix(rnorm(N * K), nrow = N)
    Z %*% L + mu_mat
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
      seq8    <- alternating_blocks()
      order10 <- c("H","L", seq8)
      env_prob[[i]] <- lapply(order10, function(tag) if (tag=="H") c(1,0) else c(0,1))
    }
  }
  
  samples_per_board <- t(replicate(
    N_id, random_session_counts(N_board, min_board, max_board, total_samples)
  ))
  
  # ---------- sampling phase counts: ALWAYS multinomial ----------
  y_sampling <- array(0L, dim = c(N_id, N_board, N_food))
  for (i in seq_len(N_id)) for (b in seq_len(N_board)) {
    size <- samples_per_board[i,b]; p <- env_prob[[i]][[b]]
    k1 <- rbinom(1, size, p[1]); y_sampling[i,b,1] <- k1; y_sampling[i,b,2] <- size-k1
  }
  
  # ---------- initial OH for test ----------
  alpha_s_vec <- ifelse(treatment == "mixed", alpha_mixed, alpha_uniform)
  
  # group-specific β (simplex)
  stopifnot(length(beta_mixed) == 2, length(beta_uniform) == 2)
  beta_mixed   <- normalize(beta_mixed)
  beta_uniform <- normalize(beta_uniform)
  
  beta_s_mat <- matrix(NA_real_, nrow = N_id, ncol = 2)
  for (i in seq_len(N_id)) {
    beta_s_mat[i, ] <- if (treatment[i] == "mixed") beta_mixed else beta_uniform
  }
  
  # ---------- per-ID random effects for (λ, φ, κ) with treatment-specific RAW means ----------
  # Convert OUTCOME-scale means (L, P, K) to RAW means (log L, logit P, log K)
  mu_by_id_raw <- matrix(NA_real_, nrow = N_id, ncol = 3)
  for (i in seq_len(N_id)) {
    if (treatment[i] == "mixed") {
      mu_by_id_raw[i, ] <- c(log(lambda_mixed),  logit(phi_mixed),  log(pmax(kappa_mixed, .Machine$double.eps)))
    } else { # uniform
      mu_by_id_raw[i, ] <- c(log(lambda_uniform), logit(phi_uniform), log(pmax(kappa_uniform, .Machine$double.eps)))
    }
  }
  colnames(mu_by_id_raw) <- c("lambda_raw_i","phi_raw_i","kappa_raw_i")
  
  
  # Covariance across (λ, φ, κ) in RAW space
  sds3 <- c(sd_lambda, sd_phi, sd_kappa)
  if (corr_type == "independent") {
    Corr3 <- diag(3)
  } else if (corr_type == "exchangeable") {
    rho <- rho_exchangeable
    Corr3 <- matrix(rho, 3, 3); diag(Corr3) <- 1
  } else {
    Corr3 <- corr_matrix_manual
  }
  Sigma3 <- diag(sds3) %*% Corr3 %*% diag(sds3)
  
  RE_ID_raw <- if (heterogeneity) rmvnorm_chol_means(mu_by_id_raw, Sigma3) else mu_by_id_raw
  colnames(RE_ID_raw) <- c("lambda_raw_i","phi_raw_i","kappa_raw_i")
  
  # Natural scales per ID (used in utilities & learning)
  L_id <- exp(RE_ID_raw[, "lambda_raw_i"])   # >0
  P_id <- invlogit(RE_ID_raw[, "phi_raw_i"]) # in (0,1)
  K_id <- exp(RE_ID_raw[, "kappa_raw_i"])    # >=0
  
  # warn if phi can't matter
  if (!is.null(fixed_A)) {
    warning("fixed_A supplied: attractions are held fixed; phi has no effect in the test phase.")
  }
  
  # ---------- trackers ----------
  P_S_tracker       <- matrix(NA_real_, nrow = N_id, ncol = max_test_trials)
  theta_b11_tracker <- array(NA_real_, dim = c(N_id, max_test_trials, N_food))
  theta_b12_tracker <- array(NA_real_, dim = c(N_id, max_test_trials, N_food))
  alpha_tracker     <- matrix(NA_real_, nrow = N_id, ncol = max_test_trials)
  beta_tracker      <- array(NA_real_, dim = c(N_id, max_test_trials, N_food))
  AF_tracker        <- array(NA_real_, dim = c(N_id, max_test_trials, N_food))
  
  # ---------- outputs (Stan data fields) ----------
  Switch      <- integer(0)  # 1=stay, 2=switch
  id_test     <- integer(0)
  trial_num   <- integer(0)
  y_test      <- integer(0)  # LV on B11
  y_test_ate  <- integer(0)  # 0/1
  T_id        <- integer(N_id)
  
  # ---------- test phase (exponential decay only) ----------
  for (i in seq_len(N_id)) {
    a_first <- alpha_s_vec[i]
    b_first <- as.numeric(beta_s_mat[i, ])
    
    alpha_t <- a_first
    beta_t  <- b_first
    counts_b11 <- c(0.0, 0.0)
    
    # attractions: fixed vs learning
    if (is.null(fixed_A)) {
      A_F <- as.numeric(A0)
      do_learn_A <- TRUE
    } else {
      stopifnot(length(fixed_A) == 2, all(is.finite(fixed_A)), all(fixed_A >= 0), all(fixed_A <= 1))
      A_F <- as.numeric(fixed_A)
      do_learn_A <- FALSE
    }
    
    t_done <- 0
    
    for (t in 1:max_test_trials) {
      # pseudo-counts and predictive per board
      OH <- as.numeric(alpha_t) * beta_t
      theta_b11 <- OH + counts_b11
      theta_b12 <- OH
      P_F_b11 <- normalize(theta_b11)
      P_F_b12 <- normalize(theta_b12)
      
      # utilities → stay/switch probabilities (per-ID L and K)
      U_b11 <- sum(P_F_b11 * A_F)
      U_b12 <- sum(P_F_b12 * A_F)
      util  <- c(L_id[i] * U_b11,
                 L_id[i] * U_b12 - K_id[i])
      P     <- softmax(util)  # c(P_stay, P_switch)
      
      # record a test observation BEFORE outcome
      y_test     <- c(y_test, 2L)
      ate        <- rbinom(1, 1, p_eat_low)
      y_test_ate <- c(y_test_ate, ate)
      
      # trackers (snapshot)
      theta_b11_tracker[i, t, ] <- P_F_b11
      theta_b12_tracker[i, t, ] <- P_F_b12
      alpha_tracker[i, t]       <- alpha_t
      beta_tracker[i, t, ]      <- beta_t
      AF_tracker[i, t, ]        <- A_F
      P_S_tracker[i, t]         <- P[2]
      
      # realised choice
      choice    <- sample(1:2, 1, prob = P)
      Switch    <- c(Switch, choice)
      id_test   <- c(id_test, i)
      trial_num <- c(trial_num, t)
      
      if (break_on_switch && choice == 2) { t_done <- t; break }
      
      # stayed → observe LV on b11, update counts
      counts_b11[2] <- counts_b11[2] + 1.0
      
      # attraction updating for LV (φ matters only if do_learn_A == TRUE)
      if (do_learn_A) {
        P_upd <- P_id[i]   # inv_logit(phi_raw_i)
        if (ate == 1) {
          A_F[2] <- (1 - P_upd) * A_F[2] + P_upd * 0.25
        } else {
          A_F[2] <- (1 - P_upd) * A_F[2]
        }
      }
      
      # advance α and β to next step (DECAY mode)
      alpha_t <- (1 - d_alpha) * alpha_t
      beta_t1 <- (1 - d_beta)  * beta_t[1]
      beta_t  <- c(beta_t1, 1 - beta_t1)
    }
    
    if (break_on_switch) {
      if (t_done == 0) t_done <- max_test_trials
      T_id[i] <- t_done
    } else {
      T_id[i] <- max_test_trials
    }
  }
  
  T_max <- max(T_id)
  
  # ---------- Preference test & totals ----------
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
    RE_ID_raw = RE_ID_raw,      # per-ID draws in RAW space (λ, φ, κ)
    L_id      = L_id,           # natural
    P_update  = P_id,           # natural
    K_id      = K_id,           # natural
    trackers = list(
      P_S = P_S_tracker,
      theta_b11 = theta_b11_tracker,
      theta_b12 = theta_b12_tracker,
      alpha_t = alpha_tracker,
      beta_t  = beta_tracker,
      AF_t    = AF_tracker
    ),
    # params used
    args = list(
      total_samples = total_samples,
      min_board = min_board,
      max_board = max_board,
      alpha_mixed = alpha_mixed,
      alpha_uniform = alpha_uniform,
      beta_mixed = beta_mixed,
      beta_uniform = beta_uniform,
      
      # outcome-scale means (echo)
      lambda_mixed   = lambda_mixed,
      lambda_uniform = lambda_uniform,
      phi_mixed      = phi_mixed,
      phi_uniform    = phi_uniform,
      kappa_mixed    = kappa_mixed,
      kappa_uniform  = kappa_uniform,
      
      # (optional) derived RAW means for reference
      lambda_raw_mixed   = log(lambda_mixed),
      lambda_raw_uniform = log(lambda_uniform),
      phi_raw_mixed      = log(phi_mixed/(1-phi_mixed)),
      phi_raw_uniform    = log(phi_uniform/(1-phi_uniform)),
      kappa_raw_mixed    = log(pmax(kappa_mixed, .Machine$double.eps)),
      kappa_raw_uniform  = log(pmax(kappa_uniform, .Machine$double.eps)),
      
      heterogeneity = heterogeneity,
      sd_lambda = sd_lambda,
      sd_phi    = sd_phi,
      sd_kappa  = sd_kappa,
      corr_type = corr_type,
      rho_exchangeable = rho_exchangeable,
      corr_matrix_manual = corr_matrix_manual,
      max_test_trials = max_test_trials,
      break_on_switch = break_on_switch,
      A0 = A0,
      fixed_A = fixed_A,
      p_eat_low = p_eat_low,
      p_eat_high = p_eat_high,
      d_alpha = d_alpha,
      d_beta  = d_beta
    )
  )
}
