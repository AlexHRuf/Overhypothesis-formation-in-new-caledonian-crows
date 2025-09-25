#' Tidy extraction of key parameters and time-varying outputs from the combo Stan fit,
#' but ONLY for parameters that are present in the fit.
#'
#' @param fit      rstan::stanfit (or compatible with tidybayes::spread_draws)
#' @param dataList the exact list you passed to Stan (used for id/trial lookups)
#' @param width    interval width for mean_qi summaries (default 0.80)
#' @return list with raw tidy draws (data frames that exist), summaries, and meta
extract_OH_draws <- function(fit, dataList, width = 0.80) {
  # --- required pkgs (use namespaces; do not attach) ---
  if (!requireNamespace("tidybayes", quietly = TRUE)) stop("Please install {tidybayes}.")
  if (!requireNamespace("dplyr",     quietly = TRUE)) stop("Please install {dplyr}.")
  if (!requireNamespace("tibble",    quietly = TRUE)) stop("Please install {tibble}.")
  if (!requireNamespace("rstan",     quietly = TRUE)) stop("Please install {rstan}.")
  if (!requireNamespace("rlang",     quietly = TRUE)) stop("Please install {rlang}.")
  if (!requireNamespace("tidyselect", quietly = TRUE)) stop("Please install {tidyselect}.")
  
  
  dplyr  <- asNamespace("dplyr")
  tb     <- asNamespace("tidybayes")
  tibble <- asNamespace("tibble")
  rlang  <- asNamespace("rlang")
  tidyselect  <- asNamespace("tidyselect")
  
  
  # --- discover parameter names in fit ---
  arr <- rstan::extract(fit, permuted = FALSE, inc_warmup = FALSE)
  param_names <- dimnames(arr)[[3]]
  
  present <- function(name) any(grepl(paste0("^", name, "(\\[|$)"), param_names))
  maybe_spread <- function(expr_present, spread_call) {
    if (expr_present) eval(spread_call) else NULL
  }
  
  # -------- RAW DRAWS (conditional) --------
  scalar_params <- c("PT_alpha","PT_weight","ALL_alpha","ALL_weight",
                     "lambda_raw","phi_raw","d_alpha","d_beta","kappa_raw")
  scalar_params <- scalar_params[vapply(scalar_params, present, logical(1))]
  scalars_df <- if (length(scalar_params)) {
    tb$spread_draws(fit, !!!rlang$syms(scalar_params))
  } else NULL
  
  alpha_df <- maybe_spread(present("alpha_s"),
                           quote(tb$spread_draws(fit, alpha_s[id])))
  
  beta_df  <- maybe_spread(present("beta_s"),
                           quote(tb$spread_draws(fit, beta_s[id, food])))
  
  PT_F_df  <- maybe_spread(present("PT_F_value"),
                           quote(tb$spread_draws(fit, PT_F_value[id, food])))
  
  ALL_F_df <- maybe_spread(present("ALL_F_value"),
                           quote(tb$spread_draws(fit, ALL_F_value[id, food])))
  
  probS_df_raw <- maybe_spread(present("Prob_S"),
                               quote(tb$spread_draws(fit, Prob_S[test, choice])))
  
  a_t_df  <- maybe_spread(present("a_t"),
                          quote(tb$spread_draws(fit, a_t[id, trial])))
  
  b_t_df  <- maybe_spread(present("b_t"),
                          quote(tb$spread_draws(fit, b_t[id, food, trial])))
  
  AF_t_df <- maybe_spread(present("AF_t"),
                          quote(tb$spread_draws(fit, AF_t[id, food, trial])))
  
  # -------- LABELS / FACTORS (only if those frames exist) --------
  food_levels <- if (!is.null(dataList$food_names)) dataList$food_names else c("high value", "low value")
  lvl_idx     <- seq_along(food_levels)
  
  if (!is.null(beta_df))  beta_df  <- dplyr$mutate(beta_df,  food = factor(.data$food, levels = lvl_idx, labels = food_levels))
  if (!is.null(PT_F_df))  PT_F_df  <- dplyr$mutate(PT_F_df,  food = factor(.data$food, levels = lvl_idx, labels = food_levels))
  if (!is.null(ALL_F_df)) ALL_F_df <- dplyr$mutate(ALL_F_df, food = factor(.data$food, levels = lvl_idx, labels = food_levels))
  if (!is.null(b_t_df))   b_t_df   <- dplyr$mutate(b_t_df,   food = factor(.data$food, levels = lvl_idx, labels = food_levels))
  if (!is.null(AF_t_df))  AF_t_df  <- dplyr$mutate(AF_t_df,  food = factor(.data$food, levels = lvl_idx, labels = food_levels))
  
  # -------- TEST ROW → id & trial mapping (only if Prob_S present) --------
  probS_df <- NULL
  if (!is.null(probS_df_raw)) {
    if (is.null(dataList$id) || length(dataList$id) != length(dataList$y_test)) {
      stop("dataList must include `id` (length N_test) to map Prob_S rows to IDs.")
    }
    test_id <- as.integer(dataList$id)
    trial_per_row <- ave(test_id, test_id, FUN = seq_along)
    
    treatment_vec <- if (!is.null(dataList$treatment)) dataList$treatment else rep("all", dataList$N_id)
    
    probS_df <- dplyr$filter(probS_df_raw, .data$choice == 2L)
    probS_df <- dplyr$mutate(
      probS_df,
      id        = test_id[.data$test],
      trial     = trial_per_row[.data$test],
      treatment = factor(treatment_vec[.data$id])
    )
    probS_df <- dplyr$transmute(probS_df, .data$.draw, .data$id, .data$trial, .data$treatment, ps = .data$Prob_S)
  }
  
  # -------- PER-ID METADATA --------
  id_key <- tibble$tibble(
    id        = seq_len(dataList$N_id),
    T_id      = as.integer(dataList$T_id),
    treatment = factor(if (!is.null(dataList$treatment)) dataList$treatment else rep("all", dataList$N_id))
  )
  
  # -------- NEW: reconstruct per-ID RAW params (lambda_raw_id, phi_raw_id, kappa_raw_id if present) --------
  re_id_raw_df <- NULL
  K <- NA_integer_
  if (present("z_ID") && present("sigma_ID") && present("Rho_ID")) {
    # Grab draws for RE pieces
    z_df     <- tb$spread_draws(fit, z_ID[k, id])
    sigma_df <- tb$spread_draws(fit, sigma_ID[k])
    R_df     <- tb$spread_draws(fit, Rho_ID[r, c])
    
    # population means (scalars)
    lam_mean <- if (present("lambda_raw")) tb$spread_draws(fit, lambda_raw) else NULL
    phi_mean <- if (present("phi_raw"))    tb$spread_draws(fit, phi_raw)    else NULL
    kap_mean <- if (present("kappa_raw"))  tb$spread_draws(fit, kappa_raw)  else NULL
    
    K <- max(z_df$k)
    N <- dataList$N_id
    
    # Build per-draw matrices
    re_id_raw_df <- z_df |>
      dplyr$group_by(.data$.draw) |>
      dplyr$group_modify(function(dd, key) {
        # Z: K x N (ordered by k, id)
        Z <- matrix(dd$z_ID, nrow = K, ncol = N, byrow = FALSE)
        # sigma: length K
        s <- sigma_df$sigma_ID[sigma_df$.draw == key$.draw]
        # Rho lower Cholesky: K x K
        Rsub <- R_df[R_df$.draw == key$.draw, ]
        Lcorr <- matrix(0, K, K)
        Lcorr[cbind(Rsub$r, Rsub$c)] <- Rsub$Rho_ID  # r,c are 1-based in tidybayes
        # M = diag(sigma) %*% Lcorr
        M <- diag(s) %*% Lcorr
        V <- t(M %*% Z)  # N x K ; rows = id, cols = dims
        
        out <- tibble::tibble(
          id = seq_len(N),
          v1 = V[, 1],
          v2 = if (K >= 2) V[, 2] else NA_real_,
          v3 = if (K >= 3) V[, 3] else NA_real_
        )
        out
      }) |>
      dplyr$ungroup()
    
    # Add population means and compute per-ID raw params
    if (!is.null(lam_mean)) re_id_raw_df <- dplyr$left_join(re_id_raw_df, lam_mean, by = ".draw")
    if (!is.null(phi_mean)) re_id_raw_df <- dplyr$left_join(re_id_raw_df, phi_mean, by = ".draw")
    if (!is.null(kap_mean)) re_id_raw_df <- dplyr$left_join(re_id_raw_df, kap_mean, by = ".draw")
    
    # ---- NEW: ensure columns exist; if missing, fill with 0 so "+ v*" works
    if (!("lambda_raw" %in% names(re_id_raw_df))) re_id_raw_df$lambda_raw <- 0
    if (!("phi_raw"    %in% names(re_id_raw_df))) re_id_raw_df$phi_raw    <- 0
    if (!("kappa_raw"  %in% names(re_id_raw_df))) re_id_raw_df$kappa_raw  <- 0
    
    # Compute per-ID raw params, then attach treatment
    re_id_raw_df <- re_id_raw_df |>
      dplyr$mutate(
        lambda_raw_id = .data$lambda_raw + .data$v1,
        phi_raw_id    = .data$phi_raw    + .data$v2,
        kappa_raw_id  = if (K >= 3) .data$kappa_raw + .data$v3 else NA_real_
      ) |>
      dplyr$select(.data$.draw, .data$id, .data$lambda_raw_id, .data$phi_raw_id, .data$kappa_raw_id) |>
      dplyr$left_join(id_key, by = "id")
  }
  # -------- SUMMARIES (only for existing pieces) --------
  summaries <- list()
  
  if (!is.null(alpha_df)) {
    sum_alpha <- dplyr$group_by(alpha_df, .data$id)
    sum_alpha <- dplyr$summarise(
      sum_alpha,
      mean = mean(.data$alpha_s),
      qlo  = stats::quantile(.data$alpha_s, (1 - width)/2),
      qhi  = stats::quantile(.data$alpha_s, 1 - (1 - width)/2),
      .groups = "drop"
    )
    sum_alpha <- dplyr$left_join(sum_alpha, id_key, by = "id")
    summaries$sum_alpha <- sum_alpha
  }
  
  if (!is.null(beta_df)) {
    sum_beta <- dplyr$group_by(beta_df, .data$id, .data$food)
    sum_beta <- dplyr$summarise(
      sum_beta,
      mean = mean(.data$beta_s),
      qlo  = stats::quantile(.data$beta_s, (1 - width)/2),
      qhi  = stats::quantile(.data$beta_s, 1 - (1 - width)/2),
      .groups = "drop"
    )
    sum_beta <- dplyr$left_join(sum_beta, id_key, by = "id")
    summaries$sum_beta <- sum_beta
  }
  
  if (!is.null(ALL_F_df)) {
    tmp <- dplyr$mutate(ALL_F_df, attr = plogis(.data$ALL_F_value))
    sum_ALL_F <- dplyr$group_by(tmp, .data$id, .data$food)
    sum_ALL_F <- dplyr$summarise(
      sum_ALL_F,
      mean = mean(.data$attr),
      qlo  = stats::quantile(.data$attr, (1 - width)/2),
      qhi  = stats::quantile(.data$attr, 1 - (1 - width)/2),
      .groups = "drop"
    )
    sum_ALL_F <- dplyr$left_join(sum_ALL_F, id_key, by = "id")
    summaries$sum_ALL_F <- sum_ALL_F
  }
  
  if (!is.null(probS_df)) {
    sum_probS <- tidybayes::mean_qi(dplyr::group_by(probS_df, id, trial, treatment), ps, .width = width)
    sum_probS <- dplyr::rename(sum_probS, mean = ps)
    summaries$sum_probS <- sum_probS
  }
  
  if (!is.null(a_t_df)) {
    tmp <- dplyr$left_join(a_t_df, id_key, by = "id")
    tmp <- dplyr$filter(tmp, .data$trial <= .data$T_id)
    tmp <- dplyr$group_by(tmp, .data$id, .data$trial, .data$treatment)
    sum_a_t <- tidybayes::mean_qi(tmp, a_t, .width = width)
    sum_a_t <- dplyr::rename(sum_a_t, mean = a_t)
    summaries$sum_a_t <- sum_a_t
  }
  
  if (!is.null(b_t_df)) {
    tmp <- dplyr$left_join(b_t_df, id_key, by = "id")
    tmp <- dplyr$filter(tmp, .data$trial <= .data$T_id)
    tmp <- dplyr$group_by(tmp, .data$id, .data$trial, .data$food, .data$treatment)
    sum_b_t <- tidybayes::mean_qi(tmp, b_t, .width = width)
    sum_b_t <- dplyr::rename(sum_b_t, mean = b_t)
    summaries$sum_b_t <- sum_b_t
  }
  
  if (!is.null(AF_t_df)) {
    tmp <- dplyr$left_join(AF_t_df, id_key, by = "id")
    tmp <- dplyr$filter(tmp, .data$trial <= .data$T_id)
    tmp <- dplyr$group_by(tmp, .data$id, .data$trial, .data$food, .data$treatment)
    sum_AF_t <- tidybayes::mean_qi(tmp, AF_t, .width = width)
    sum_AF_t <- dplyr$rename(sum_AF_t, mean = AF_t)
    summaries$sum_AF_t <- sum_AF_t
  }
  
  # NEW: summaries for per-ID raws
  if (!is.null(re_id_raw_df)) {
    have_kappa <- "kappa_raw_id" %in% names(re_id_raw_df)
    cols <- c("lambda_raw_id", "phi_raw_id", if (have_kappa) "kappa_raw_id")
    
    sum_re <- re_id_raw_df |>
      dplyr$group_by(.data$id, .data$treatment) |>
      dplyr$summarise(
        dplyr$across(
          tidyselect::all_of(cols),
          list(
            mean = ~ mean(.x, na.rm = TRUE),
            qlo  = ~ stats::quantile(.x, (1 - width)/2, na.rm = TRUE),
            qhi  = ~ stats::quantile(.x, 1 - (1 - width)/2, na.rm = TRUE)
          ),
          .names = "{.col}_{.fn}"
        ),
        .groups = "drop"
      )
    
    summaries$re_id_raw <- sum_re
  }
  
  # assemble draws list
  draws <- list()
  if (!is.null(scalars_df))   draws$scalars   <- scalars_df
  if (!is.null(alpha_df))     draws$alpha     <- alpha_df
  if (!is.null(beta_df))      draws$beta      <- beta_df
  if (!is.null(PT_F_df))      draws$PT_F      <- PT_F_df
  if (!is.null(ALL_F_df))     draws$ALL_F     <- ALL_F_df
  if (!is.null(probS_df))     draws$probS     <- probS_df
  if (!is.null(a_t_df))       draws$a_t       <- a_t_df
  if (!is.null(b_t_df))       draws$b_t       <- b_t_df
  if (!is.null(AF_t_df))      draws$AF_t      <- AF_t_df
  if (!is.null(re_id_raw_df)) draws$re_id_raw <- re_id_raw_df
  
  # -------- RETURN --------
  list(
    meta = list(
      param_names = param_names,
      food_levels = food_levels,
      id_key      = id_key,
      width       = width,
      RE_dim      = K
    ),
    draws = draws,
    summaries = summaries
  )
}
