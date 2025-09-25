traceplot_first <- function(fit,
                            n_keep = 500,
                            color_scheme = "brewer-Spectral",
                            include_lp = TRUE) {
  # deps
  if (!requireNamespace("bayesplot", quietly = TRUE))
    stop("Package 'bayesplot' is required.")
  if (!requireNamespace("rstan", quietly = TRUE))
    stop("Package 'rstan' is required.")
  
  # set colors (won't error if unknown; bayesplot has several brewer-* schemes)
  try(bayesplot::color_scheme_set(color_scheme), silent = TRUE)
  
  # 1) Full parameter names with indices (e.g., "alpha_s[1]", "beta_s[1,1]")
  arr <- as.array(fit)  # chains x iterations x parameters
  dn  <- dimnames(arr)
  if (is.null(dn) || length(dn) < 3 || is.null(dn[[3]]))
    stop("Couldn't retrieve parameter names via dimnames(as.array(fit)).")
  full_names <- dn[[3]]
  
  # 2) Base families from rstan::extract (list names, e.g., "alpha_s", "beta_s")
  fams <- names(rstan::extract(fit))
  if (is.null(fams)) fams <- character(0)
  
  # helper: first occurrence of a family's indexed parameter
  pick_first <- function(par) {
    if (par %in% full_names) return(par) # scalar parameter (no indices)
    # starts with "par[" (use grepl to be safe)
    hits <- full_names[grepl(paste0("^", gsub("\\[", "\\\\[", par), "\\["), full_names)]
    if (length(hits)) hits[1] else NA_character_
  }
  
  pars_first <- unique(stats::na.omit(vapply(fams, pick_first, character(1))))
  if (include_lp && "lp__" %in% full_names) {
    pars_first <- unique(c(pars_first, "lp__"))
  }
  
  if (length(pars_first) == 0L)
    stop("No parameters found to plot. Check that 'fit' has named parameters.")
  
  # 3) Thinning for plotting
  iter <- dim(arr)[2]
  keep <- sort(unique(round(seq(1, iter, length.out = min(n_keep, iter)))))
  arr_small <- arr[, keep, , drop = FALSE]
  
  # 4) Plot
  bayesplot::mcmc_trace(arr_small, pars = pars_first)
}
