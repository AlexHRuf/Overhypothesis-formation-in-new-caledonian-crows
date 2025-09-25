generate_stan_data_simple <- function(
    # ---- sizes ----
    N_id = 16,
    N_board = 10,
    N_food = 2,
    seed = 12345,
    
    # ---- group assignment ----
    treatment = NULL,  # if NULL, auto-generate alternating "mixed"/"uniform"
    
    # ---- sampling environment ----
    total_samples = 120,
    min_board     = 8,
    max_board     = 14,
    evidence_noise = c("multinomial","deterministic")[1],
    
    # ---- per-group test-length controls ----
    mean_T_mixed   = 50,
    mean_T_uniform = 45,
    sd_T           = 5,
    T_min          = 1,
    T_max_cap      = 50,
    
    # ---- probabilities you set ----
    p_mixed_board = c(0.5, 0.5),
    p_eat_high = 0.95,
    p_eat_low  = 0.05,
    
    # preference test meta
    PT_trials_range  = c(60L, 80L),
    PT_log_dur_range = c(log(300), log(600))
) {
  stopifnot(N_food == 2)
  set.seed(seed)
  
  # ---------- handle treatment ----------
  if (is.null(treatment)) {
    treatment <- rep(c("mixed","uniform"), length.out = N_id)
  } else {
    if (length(treatment) == 1) {
      treatment <- rep(treatment, N_id)
    } else if (length(treatment) != N_id) {
      stop("treatment must be length 1, N_id, or NULL.")
    }
  }
  # ---------- helpers ----------
  normalize <- function(v) { v <- as.numeric(v); s <- sum(v); if (!is.finite(s) || s <= 0) rep(1/length(v), length(v)) else v/s }
  
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
    # length-8 sequence of H/L with no >2 identical in a row
    v <- c(rep("H",4), rep("L",4))
    tries <- 0
    repeat {
      w <- sample(v, length(v), replace = FALSE)
      if (all(rle(w)$lengths <= 2)) return(w)
      tries <- tries + 1
      if (tries > 5000) stop("Failed to make uniform sequence with <=2 in a row.")
    }
  }
  
  # ---------- board-level env probs by treatment ----------
  env_prob <- vector("list", N_id)
  for (i in seq_len(N_id)) {
    if (treatment[i] == "mixed") {
      env_prob[[i]] <- replicate(N_board, normalize(p_mixed_board), simplify = FALSE)
    } else { # uniform: pure H or L boards with no long runs
      seq8    <- alternating_blocks()
      order10 <- c("H","L", seq8)  # length 10
      env_prob[[i]] <- lapply(order10, function(tag) if (tag=="H") c(1,0) else c(0,1))
    }
  }
  
  # ---------- per-ID per-board sample sizes ----------
  samples_per_board <- t(replicate(
    N_id, random_session_counts(N_board, min_board, max_board, total_samples)
  ))
  
  # ---------- y_sampling generation ----------
  y_sampling <- array(0L, dim = c(N_id, N_board, N_food))
  for (i in seq_len(N_id)) for (b in seq_len(N_board)) {
    size <- samples_per_board[i, b]
    p    <- env_prob[[i]][[b]]
    if (evidence_noise == "multinomial") {
      k1 <- rbinom(1, size, p[1])
      y_sampling[i, b, 1] <- k1
      y_sampling[i, b, 2] <- size - k1
    } else { # deterministic rounding
      y_sampling[i, b, ] <- as.integer(round(p * size))
      diff <- size - sum(y_sampling[i, b, ])
      if (diff != 0) { j <- which.max(p); y_sampling[i, b, j] <- y_sampling[i, b, j] + diff }
    }
  }
  
  # ---------- automatic per-group T_id drawing ----------
  mu_by_group <- ifelse(treatment == "mixed", mean_T_mixed, mean_T_uniform)
  T_id <- round(rnorm(N_id, mean = mu_by_group, sd = sd_T))
  T_id <- pmin(pmax(T_id, T_min), T_max_cap)
  T_id <- as.integer(T_id)
  T_max <- as.integer(max(T_id))
  
  # ---------- preference test totals from your probs ----------
  PT_trials  <- sample(PT_trials_range[1]:PT_trials_range[2], N_id, TRUE)
  PT_log_dur <- runif(N_id, PT_log_dur_range[1], PT_log_dur_range[2])
  
  PT_ate_HV  <- rbinom(n = N_id, size = PT_trials, prob = p_eat_high)
  PT_ate_LV  <- rbinom(n = N_id, size = PT_trials, prob = p_eat_low)
  PT_total_eat <- cbind(PT_ate_HV, PT_ate_LV)
  
  # ---------- ALL_total_eat from sampling exposures + your probs ----------
  samp_counts <- apply(y_sampling, c(1,3), sum) # N_id x 2
  samp_ate_HV <- rbinom(n = N_id, size = samp_counts[,1], prob = p_eat_high)
  samp_ate_LV <- rbinom(n = N_id, size = samp_counts[,2], prob = p_eat_low)
  samp_ate <- cbind(samp_ate_HV, samp_ate_LV)
  
  ALL_total_eat <- PT_total_eat + samp_ate
  ALL_log_dur   <- PT_log_dur * 2
  
  # ---------- TEST PHASE (always auto) ----------
  # Build per-ID sequences: y_test=2; Switch=1 except last trial -> 2; y_test_ate ~ Bernoulli(p_eat_low)
  N_test <- sum(T_id)
  id_vec <- inverse.rle(list(values = seq_len(N_id), lengths = T_id))
  
  y_test <- rep(2L, N_test)
  y_test_ate <- rbinom(N_test, 1, p_eat_low)
  
  Switch <- integer(N_test)
  idx <- 1L
  for (j in seq_len(N_id)) {
    tlen <- T_id[j]
    if (tlen <= 0L) next
    if (tlen == 1L) {
      Switch[idx] <- 2L
      idx <- idx + 1L
    } else {
      Switch[idx:(idx + tlen - 2L)] <- 1L
      Switch[idx + tlen - 1L] <- 2L
      idx <- idx + tlen
    }
  }
  
  # ---------- coerce to Stan-friendly types ----------
  as_int <- function(x) { storage.mode(x) <- "integer"; x }
  
  list(
    # sizes / indexing
    N_test  = as_int(N_test),
    N_food  = as_int(N_food),
    N_id    = as_int(N_id),
    N_board = as_int(N_board),
    T_id    = as_int(T_id),
    T_max   = as_int(T_max),
    
    # preference test (exposure)
    PT_total_eat = as_int(PT_total_eat),
    PT_log_dur   = PT_log_dur,
    
    # sampling phase
    y_sampling = as_int(y_sampling),
    
    # pref + sampling totals (exposure)
    ALL_total_eat = as_int(ALL_total_eat),
    ALL_log_dur   = ALL_log_dur,
    
    # test phase (auto)
    y_test      = as_int(y_test),
    y_test_ate  = as_int(y_test_ate),
    id          = as_int(id_vec),
    Switch      = as_int(Switch),
    
    # ------ extras (not used by Stan, but handy) ------
    treatment   = treatment,
    samples_per_board = samples_per_board,
    args = list(
      seed = seed,
      total_samples = total_samples,
      min_board = min_board,
      max_board = max_board,
      evidence_noise = evidence_noise,
      mean_T_mixed = mean_T_mixed,
      mean_T_uniform = mean_T_uniform,
      sd_T = sd_T,
      T_min = T_min,
      T_max_cap = T_max_cap,
      p_mixed_board = p_mixed_board,
      p_eat_high = p_eat_high,
      p_eat_low  = p_eat_low,
      PT_trials_range = PT_trials_range,
      PT_log_dur_range = PT_log_dur_range
    )
  )
}
