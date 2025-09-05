#' .calc_mic (internal)
#' Compute LP‑scale MIC, transformed MIC, gradients, and SEs.
#' @keywords internal
#' @noRd
.calc_mic <- function(all_coefs, vcov_mat, alpha1_name, beta_c_terms,
                      beta_nc_terms, mm_c, mm_nc, inv_fun) {
  alpha1_val <- all_coefs[alpha1_name]
  numer      <- as.numeric(alpha1_val - mm_nc %*% all_coefs[beta_nc_terms])
  denom      <- as.numeric(mm_c       %*% all_coefs[beta_c_terms])
  lp         <- numer / denom
  mic        <- inv_fun(lp)

  g_lp <- matrix(0, nrow = length(lp), ncol = length(all_coefs),
                 dimnames = list(NULL, names(all_coefs)))
  g_lp[, alpha1_name] <- 1 / denom
  if (length(beta_nc_terms)) g_lp[, beta_nc_terms] <- -mm_nc / denom
  g_lp[, beta_c_terms] <- mm_c * (-lp / denom)

  d_mic_d_lp <- if (identical(inv_fun, expm1)) (mic + 1) else exp(lp)
  g_mic      <- g_lp * d_mic_d_lp
  # se_lp      <- sqrt(pmax(0, rowSums((g_lp %*% vcov_mat) * g_lp)))
  g_lpA <- .align_grad_mat_to_vcov(g_lp, vcov_mat)         # n×k, name-aligned
  se_lp <- sqrt(pmax(0, rowSums((g_lpA %*% vcov_mat) * g_lpA)))


  list(lp = lp, mic = mic, g_lp = g_lp, g_mic = g_mic, se_lp = se_lp)
}
