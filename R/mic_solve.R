#' Estimate MICs and Compare Groups
#'
#' High‑level wrapper that returns:
#' * `delta_mic_results` – additive pairwise differences (Δ‑MIC).
#' * `ratio_mic_results` – multiplicative pairwise ratios.
#' * `dod_ratio_results` – **difference‑of‑differences on the ratio scale**
#'   (ratio‑of‑ratios, a classic interaction on the log scale).
#' * `dod_delta_results` – **difference‑of‑differences on the additive scale**
#'   (Δ of Δs).
#'
#' @param clm_fit  Fitted object from [ordinal::clm()].
#' @param newdata  Data frame with factor combinations to evaluate.
#' @param conc_name Character string giving the **raw** concentration column.
#' @param transform_fun Transformation used in the model (default `log1p`).
#' @param inv_transform_fun Inverse transformation (default `expm1`).
#' @param alpha     Confidence‑level significance (default 0.05).
#' @param compare_pairs  One of \code{"all"} (default) to retain every pairwise
#'   comparison, or \code{"share_any"} to exclude contrasts where the two groups
#'   share no covariate levels in \code{newdata}.
#' @param pvalue_scale Which pivot the main P_value uses for pairwise tests:
#'   "lp" (difference in lp = log1p(MIC), recommended for calibration) or
#'   "logmic" (current Wald on log(MIC) for ratios and MIC scale for deltas).
#' @return An object of class `"mic_solve"` containing the tables above plus
#'   `mic_estimates`.
#' @examples
#' if (requireNamespace("ordinal", quietly = TRUE)) {
#'   ## Toy ordinal dataset
#'   set.seed(1)
#'   fit <- ordinal::clm(score ~ strain * treatment + log1p(conc), data = yeast_df)
#'   res <- mic_solve(fit, conc_name = "conc")
#'   head(res$ratio_mic_results)
#' }
#' @export
mic_solve <- function(clm_fit, newdata = NULL, conc_name,
                      transform_fun = log1p, inv_transform_fun = expm1,
                      alpha = 0.05, compare_pairs = "all",
                      pvalue_scale = c("lp", "logmic")) {

  pvalue_scale <- match.arg(pvalue_scale)
  ## --------------------------------------------------------------- ##
  ## if user omitted newdata, build one from the model’s factor levels
  if (is.null(newdata)) {
    if (length(clm_fit$xlevels) == 0)
      stop("'newdata' is missing and the model has no factor predictors.")
    newdata <- expand.grid(clm_fit$xlevels, stringsAsFactors = FALSE)
  }
  ## --------------------------------------------------------------- ##

  .validate_inputs(clm_fit, newdata, conc_name)

  if (identical(transform_fun, identity)) {
    conc_term <- conc_name              # look for plain "conc"
  } else {
    conc_term <- paste0(deparse(substitute(transform_fun)),
                        "(", conc_name, ")")
  }

  clm_call   <- clm_fit$call

  compare_pairs <- match.arg(compare_pairs, choices = c("all", "share_any"))

  zcrit      <- stats::qnorm(1 - alpha / 2)
  coefs      <- stats::coef(clm_fit)
  vcv        <- stats::vcov(clm_fit)
  beta_names <- names(clm_fit$beta)

  trans_name <- paste0(deparse(substitute(transform_fun)), "(", conc_name, ")")
  raw_name   <- conc_name

  if (trans_name %in% beta_names) {
    conc_term <- trans_name
  } else if (raw_name %in% beta_names) {
    conc_term <- raw_name
  } else {
    stop(sprintf(
      "Concentration term not found: looked for '%s' or '%s' in coefficients.",
      trans_name, raw_name))
  }
  is_c      <- grepl(conc_term, beta_names, fixed = TRUE)
  beta_c    <- beta_names[ is_c ]
  beta_nc   <- beta_names[!is_c ]

  mm <- .build_model_matrix(clm_fit, newdata, conc_name, inv_transform_fun(1))
  res <- .calc_mic(coefs, vcv, names(clm_fit$alpha)[1], beta_c, beta_nc,
                   mm[, beta_c,  drop = FALSE], mm[, beta_nc, drop = FALSE],
                   inv_transform_fun)

  mic_df <- cbind(newdata, MIC = res$mic, SE_LP = res$se_lp,
                  CI_Lower = inv_transform_fun(res$lp - zcrit * res$se_lp),
                  CI_Upper = inv_transform_fun(res$lp + zcrit * res$se_lp))

  # If there is fewer than 2 grouping factors, 'share_any' cannot produce pairs.
  if (identical(compare_pairs, "share_any") && length(names(newdata)) < 2L) {
    warning("compare_pairs = 'share_any' requires at least 2 grouping factors; using 'all' instead.", call. = FALSE)
    compare_pairs <- "all"
  }

  if (compare_pairs == "share_any") {
    cmb <- utils::combn(seq_len(nrow(newdata)), 2)
    share <- apply(cmb, 2, function(idx)
      any(newdata[idx[1], ] == newdata[idx[2], ]))
  } else {
    share <- NULL   # keep all
  }

  pw <- .pairwise_delta_ratio(res$mic, res$g_mic, log(res$mic),
                              res$g_mic / res$mic,
                              res$lp, res$g_lp, vcv, zcrit,
                              newdata, keep = share,
                              pvalue_scale)

  # dod_ratio  <- .difference_of_differences_ratio(log(res$mic),
  #                                                res$g_mic / res$mic,
  #                                                vcv, zcrit, mic_df)
  #
  # dod_delta  <- .difference_of_differences_additive(res$mic, res$g_mic,
  #                                                   vcv, zcrit, mic_df)
  dod_ratio <- .difference_of_differences_ratio_all(
    log(res$mic), res$g_mic / res$mic, vcv, zcrit, mic_df
  )

  dod_delta <- .difference_of_differences_additive_all(
    res$mic, res$g_mic, vcv, zcrit, mic_df
  )


  structure(
    list(
      mic_estimates   = mic_df,
      delta_mic_results   = pw$delta,
      ratio_mic_results   = pw$ratio,
      dod_delta_results   = dod_delta,
      dod_ratio_results   = dod_ratio,
      model               = clm_fit,
      newdata             = newdata,
      transform_fun       = transform_fun,
      inv_transform_fun   = inv_transform_fun,
      conc_name           = conc_name,
      trans_name          = trans_name,
      call                = clm_call,
      data                = clm_fit$model
      ),
            class = "mic_solve")
}

