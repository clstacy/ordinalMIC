#' .difference_of_differences_ratio (internal)
#' Ratio‑of‑ratios (multiplicative DoD) for 2×2 designs.
#' @keywords internal
#' @noRd
.difference_of_differences_ratio <- function(log_mic, g_log_mic,
                                             vcov_mat, zcrit, mic_df) {
   if (length(log_mic) != 4) return(NULL)
    l_int <- (log_mic[4] - log_mic[2]) - (log_mic[3] - log_mic[1])
    g_int <- (g_log_mic[4, ] - g_log_mic[2, ]) -
      (g_log_mic[3, ] - g_log_mic[1, ])
    se_i  <- sqrt(drop(t(g_int) %*% vcov_mat %*% g_int))
    ci_i  <- l_int + c(-1, 1) * zcrit * se_i
    # add comparison label
    comparison_df <- .dod_all_base(
                        mic_tbl = mic_df,
                        mic_col = "MIC"
                      )

    cbind(
       comparison_df,
       data.frame(
         Estimate = exp(l_int), SE_logDoD = se_i,
         CI_Lower = exp(ci_i[1]), CI_Upper = exp(ci_i[2]),
         P_value = 2 * stats::pnorm(-abs(l_int / se_i))
       )
    )
    # (comparison_df,
    #       data.frame(
    #           Estimate = exp(l_int), SE_logDoD = se_i,
    #           CI_Lower = exp(ci_i[1]), CI_Upper = exp(ci_i[2]),
    #           P_value = 2 * stats::pnorm(-abs(l_int / se_i))
    #           )
    #    )
 }

#' .difference_of_differences_additive (internal)
#' Additive Δ‑of‑Δs for 2×2 designs.
#' @keywords internal
#' @noRd
.difference_of_differences_additive <- function(mic, g_mic,
                                               vcov_mat, zcrit, mic_df) {
  if (length(mic) != 4) return(NULL)
   d_int <- (mic[4] - mic[2]) - (mic[3] - mic[1])
   g_int <- (g_mic[4, ] - g_mic[2, ]) - (g_mic[3, ] - g_mic[1, ])
   se_i  <- sqrt(drop(t(g_int) %*% vcov_mat %*% g_int))
   ci_i  <- d_int + c(-1, 1) * zcrit * se_i
   # add comparison label
   comparison_df <- .dod_all_base(
     mic_tbl = mic_df,
     mic_col = "MIC"
   )
   cbind(comparison_df,
         data.frame(Estimate = d_int, SE_DoD = se_i,
              CI_Lower = ci_i[1], CI_Upper = ci_i[2],
              P_value = 2 * stats::pnorm(-abs(d_int / se_i))
              )
   )
 }

