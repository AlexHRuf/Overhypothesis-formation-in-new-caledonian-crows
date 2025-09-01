#' extract_OH_draws()
#' -------------------
#' Tidy extraction of key parameters and time-varying outputs from a Stan fit
#' of the “over-hypothesis” model, plus convenient per-ID / per-trial summaries.
#'
#' @param fit      rstan::stanfit (or an object supported by tidybayes::spread_draws)
#' @param dataList the exact list you passed to Stan (used for id/trial lookups)
#' @param width    interval width for mean_qi summaries (default 0.80)
#'
#' @return list with raw tidy draws (data frames), summaries, and a bit of meta
extract_OH_draws <- function(fit, dataList, width = 0.80) {
  # --- packages we rely on (don’t attach; use ::) ---
  if (!requireNamespace("tidybayes", quietly = TRUE)) stop("Please install {tidybayes}.")
  if (!requireNamespace("dplyr", quietly = TRUE))      stop("Please install {dplyr}.")
  if (!requireNamespace("tibble", quietly = TRUE))     stop("Please install {tibble}.")
  if (!requireNamespace("rstan", quietly = TRUE))      stop("Please install {rstan}.")
  
  dplyr  <- asNamespace("dplyr")
  tb     <- asNamespace("tidybayes")
  tibble <- asNamespace("tibble")
  
  # --- light meta (handy for debugging) ---
  arr <- rstan::extract(fit, permuted = FALSE, inc_warmup = FALSE)
  param_names <- dimnames(arr)[[3]]
  
  # -------- RAW DRAWS --------
  scalars_df <- tb$spread_draws(
    fit,
    PT_alpha, PT_weight, ALL_alpha, ALL_weight,
    lambda_raw, phi_raw, alpha_last_raw, beta_last_raw, stickiness_raw
  )
  
  alpha_df <- tb$spread_draws(fit, alpha_s[id])                     # id
  beta_df  <- tb$spread_draws(fit, beta_s[id, food])                # id, food
  PT_F_df  <- tb$spread_draws(fit, PT_F_value[id, food])            # id, food
  ALL_F_df <- tb$spread_draws(fit, ALL_F_value[id, food])           # id, food
  
  delta_alpha_df <- tb$spread_draws(fit, delta_alpha_raw[id, trial])# id, trial
  delta_beta_df  <- tb$spread_draws(fit, delta_beta_raw[id, trial]) # id, trial
  
  probS_df_raw <- tb$spread_draws(fit, Prob_S[test, choice])        # test, choice
  a_t_df  <- tb$spread_draws(fit, a_t[id, trial])                   # id, trial
  b_t_df  <- tb$spread_draws(fit, b_t[id, food, trial])             # id, food, trial
  AF_t_df <- tb$spread_draws(fit, AF_t[id, food, trial])            # id, food, trial
  
  # -------- LABELS / FACTORS --------
  food_levels <- if (!is.null(dataList$food_names)) dataList$food_names else c("high value", "low value")
  lvl_idx     <- seq_along(food_levels)
  
  beta_df  <- dplyr$mutate(beta_df,  food = factor(food, levels = lvl_idx, labels = food_levels))
  PT_F_df  <- dplyr$mutate(PT_F_df,  food = factor(food, levels = lvl_idx, labels = food_levels))
  ALL_F_df <- dplyr$mutate(ALL_F_df, food = factor(food, levels = lvl_idx, labels = food_levels))
  b_t_df   <- dplyr$mutate(b_t_df,   food = factor(food, levels = lvl_idx, labels = food_levels))
  AF_t_df  <- dplyr$mutate(AF_t_df,  food = factor(food, levels = lvl_idx, labels = food_levels))
  
  # Map test row -> id, trial, treatment; keep only P(switch) = Prob_S[, choice==2]
  probS_df <- probS_df_raw |>
    dplyr$filter(choice == 2) |>
    dplyr$mutate(
      id        = dataList$id[test],
      trial     = dataList$trial[test],
      treatment = factor(dataList$treatment[id])
    ) |>
    dplyr$transmute(.draw, id, trial, treatment, ps = Prob_S)
  
  # -------- PER-ID METADATA & MASK PADDED TRIALS --------
  id_key <- tibble$tibble(
    id        = seq_len(dataList$N_id),
    T_id      = as.integer(dataList$T_id),
    treatment = factor(dataList$treatment)
  )
  
  # -------- SUMMARIES --------
  sum_alpha <- alpha_df |>
    dplyr$group_by(id) |>
    dplyr$summarise(
      mean = mean(alpha_s),
      qlo  = stats::quantile(alpha_s, (1 - width)/2),
      qhi  = stats::quantile(alpha_s, 1 - (1 - width)/2),
      .groups = "drop"
    ) |>
    dplyr$left_join(id_key, by = "id")
  
  sum_beta <- beta_df |>
    dplyr$group_by(id, food) |>
    dplyr$summarise(
      mean = mean(beta_s),
      qlo  = stats::quantile(beta_s, (1 - width)/2),
      qhi  = stats::quantile(beta_s, 1 - (1 - width)/2),
      .groups = "drop"
    ) |>
    dplyr$left_join(id_key, by = "id")
  
  sum_ALL_F <- ALL_F_df |>
    dplyr$mutate(attr = plogis(ALL_F_value)) |>
    dplyr$group_by(id, food) |>
    dplyr$summarise(
      mean = mean(attr),
      qlo  = stats::quantile(attr, (1 - width)/2),
      qhi  = stats::quantile(attr, 1 - (1 - width)/2),
      .groups = "drop"
    ) |>
    dplyr$left_join(id_key, by = "id")
  
  sum_probS <- probS_df |>
    dplyr$group_by(id, trial, treatment) |>
    tidybayes::mean_qi(ps, .width = width) |>
    dplyr::rename(mean = ps)
  
  sum_a_t <- a_t_df |>
    dplyr$left_join(id_key, by = "id") |>
    dplyr$filter(trial <= T_id) |>
    dplyr$group_by(id, trial, treatment) |>
    tidybayes::mean_qi(a_t, .width = width) |>
    dplyr::rename(mean = a_t)
  
  sum_b_t <- b_t_df |>
    dplyr$left_join(id_key, by = "id") |>
    dplyr$filter(trial <= T_id) |>
    dplyr$group_by(id, trial, food, treatment) |>
    tidybayes::mean_qi(b_t, .width = width) |>
    dplyr::rename(mean = b_t)
  
  sum_AF_t <- AF_t_df |>
    dplyr$left_join(id_key, by = "id") |>
    dplyr$filter(trial <= T_id) |>
    dplyr$group_by(id, trial, food, treatment) |>
    tidybayes::mean_qi(AF_t, .width = width) |>
    dplyr::rename(mean = AF_t)
  
  # -------- RETURN --------
  list(
    meta = list(
      param_names = param_names,
      food_levels = food_levels,
      id_key = id_key,
      width = width
    ),
    draws = list(
      scalars = scalars_df,
      alpha   = alpha_df,
      beta    = beta_df,
      PT_F    = PT_F_df,
      ALL_F   = ALL_F_df,
      delta_alpha = delta_alpha_df,
      delta_beta  = delta_beta_df,
      probS   = probS_df,
      a_t     = a_t_df,
      b_t     = b_t_df,
      AF_t    = AF_t_df
    ),
    summaries = list(
      sum_alpha = sum_alpha,
      sum_beta  = sum_beta,
      sum_ALL_F = sum_ALL_F,
      sum_probS = sum_probS,
      sum_a_t   = sum_a_t,
      sum_b_t   = sum_b_t,
      sum_AF_t  = sum_AF_t
    )
  )
}
