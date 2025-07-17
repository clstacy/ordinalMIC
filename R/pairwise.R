#' .pairwise_delta_ratio (internal)
#' Build Δ‑MIC and ratio tables with SEs and CIs.
#' Optionally filter out pairs that share **no** covariates in `newdata`.
#'
#' @param keep Logical vector the same length as `cmb` indicating which
#'   comparisons to keep (pre‑computed) — OR `NULL`, in which case no filtering
#'   is done.  This lightweight design avoids repeated factor comparisons.
#' @keywords internal
#' @noRd
.pairwise_delta_ratio <- function(mic, g_mic, log_mic, g_log_mic,
                                  vcov_mat, zcrit, newdata,
                                  keep = NULL) {
  cmb <- utils::combn(seq_along(mic), 2)
  if (is.null(keep)) keep <- rep(TRUE, ncol(cmb))

  out_delta <- vector("list", ncol(cmb))
  out_ratio <- vector("list", ncol(cmb))
  out_idx   <- 1L

  for (k in seq_len(ncol(cmb))) {
    if (!keep[k]) next                # skip unwanted pairs

    i <- cmb[1, k]; j <- cmb[2, k]
    g1 <- paste0(unlist(newdata[i, ]), collapse = ":")
    g2 <- paste0(unlist(newdata[j, ]), collapse = ":")

    ## Δ‑MIC
    d       <- mic[j] - mic[i]
    g_d     <- g_mic[j, ] - g_mic[i, ]
    se_d    <- sqrt(drop(t(g_d) %*% vcov_mat %*% g_d))
    ci_d    <- d + c(-1, 1) * zcrit * se_d
    p_d     <- 2 * stats::pnorm(-abs(d / se_d))
    out_delta[[k]] <- data.frame(Group1 = g1, Group2 = g2, Delta_MIC = d,
                                 SE = se_d, CI_Lower = ci_d[1],
                                 CI_Upper = ci_d[2], P_value = p_d)

    ## Ratio
    lr      <- log_mic[j] - log_mic[i]
    g_lr    <- g_log_mic[j, ] - g_log_mic[i, ]
    se_lr   <- sqrt(drop(t(g_lr) %*% vcov_mat %*% g_lr))
    ci_lr   <- lr + c(-1, 1) * zcrit * se_lr
    p_lr    <- 2 * stats::pnorm(-abs(lr / se_lr))
    out_ratio[[k]] <- data.frame(Group1 = g1, Group2 = g2,
                                 Ratio_MIC = exp(lr), SE_logRatio = se_lr,
                                 CI_Lower = exp(ci_lr[1]),
                                 CI_Upper = exp(ci_lr[2]), P_value = p_lr)

    out_idx <- out_idx + 1L
  }
  list(delta = do.call(rbind, out_delta), ratio = do.call(rbind, out_ratio))
}
