# helper to construct a single comparison label row from a 2x2 mic_df
# Assumes mic_df has exactly 4 rows and 3 columns: var1, var2, MIC (names don't matter)
# Returns data.frame with columns:
#   var1, var2, var1_lvlA, var1_lvlB, var2_lvlC, var2_lvlD
# Build a single comparison label from a 4-row slice that may contain extra constant columns.
# - Ignores columns that don't vary
# - Picks the first two columns that vary with exactly 2 levels each
# - Returns NULL if <2 eligible columns (so callers can no-op / return empty DoD)
.dod_label_from_mic_df <- function(mic_df, mic_col = "MIC") {
  if (!is.data.frame(mic_df)) stop("ordinalMIC: mic_df must be a data.frame")
  if (!mic_col %in% names(mic_df)) stop("ordinalMIC: missing MIC column: ", mic_col)

  key_cols <- setdiff(names(mic_df), mic_col)
  if (!length(key_cols)) return(NULL)

  df <- mic_df[, c(key_cols, mic_col), drop = FALSE]

  # count distinct values per design column
  nlev <- vapply(key_cols, function(k) {
    length(unique(df[[k]]))
  }, integer(1))

  # eligible: columns that vary and have exactly 2 levels
  eligible <- key_cols[nlev == 2L]

  # If more than two eligible, take the first two (consistent/stable choice)
  if (length(eligible) < 2L) {
    # Not a 2x2 slice (e.g., only one factor varies or none) -> signal no DoD
    return(NULL)
  }
  v1 <- eligible[1]; v2 <- eligible[2]

  # Coerce to factors for ordered levels
  if (!is.factor(df[[v1]])) df[[v1]] <- factor(df[[v1]])
  if (!is.factor(df[[v2]])) df[[v2]] <- factor(df[[v2]])

  # Basic sanity: the 2x2 cross should have up to 4 rows
  # (If the slice is larger, upstream should subset; we still try to label.)
  lv1 <- levels(df[[v1]])
  lv2 <- levels(df[[v2]])

  data.frame(
    var1       = v1,
    var2       = v2,
    var1_lvlA  = as.character(lv1[1]),
    var1_lvlB  = as.character(lv1[2]),
    var2_lvlC  = as.character(lv2[1]),
    var2_lvlD  = as.character(lv2[2]),
    stringsAsFactors = FALSE
  )
}


# Produce an empty (0-row) result with the SAME columns your callers expect
.dod_empty_row <- function(kind = c("ratio","additive")) {
  kind <- match.arg(kind)
  base_cols <- data.frame(
    var1=character(), var2=character(),
    var1_lvlA=character(), var1_lvlB=character(),
    var2_lvlC=character(), var2_lvlD=character(),
    stringsAsFactors = FALSE
  )
  if (kind == "ratio") {
    cbind(base_cols, data.frame(
      DDlog2MIC = numeric(), SE_logDoD = numeric(),
      CI_Lower = numeric(),  CI_Upper = numeric(),
      P_value  = numeric(),
      stringsAsFactors = FALSE
    ))
  } else {
    cbind(base_cols, data.frame(
      DDMIC = numeric(), SE_DoD = numeric(),
      CI_Lower = numeric(), CI_Upper = numeric(),
      P_value  = numeric(),
      stringsAsFactors = FALSE
    ))
  }
}

# ---------- MULTIPLICATIVE (ratio-of-ratios on log scale) ----------
#' .difference_of_differences_ratio (internal)
#' Ratio-of-ratios (multiplicative DoD) for 2x2 designs.
#' @keywords internal
#' @noRd
.difference_of_differences_ratio <- function(log_mic, g_log_mic, vcov_mat, zcrit, mic_df) {
  lab <- .dod_label_from_mic_df(mic_df, mic_col = "MIC")
  if (is.null(lab) || length(log_mic) != 4L) return(.dod_empty_row("ratio"))

  if (!is.matrix(g_log_mic) || nrow(g_log_mic) != 4L)
    stop("ordinalMIC: g_log_mic must be 4xp")
  if (!is.matrix(vcov_mat))
    stop("ordinalMIC: vcov_mat must be a square matrix")

  # (4-2) - (3-1) contrast on log scale
  l_int <- (log_mic[4] - log_mic[2]) - (log_mic[3] - log_mic[1])
  g_int <- (g_log_mic[4,,drop=FALSE] - g_log_mic[2,,drop=FALSE]) -
    (g_log_mic[3,,drop=FALSE] - g_log_mic[1,,drop=FALSE])

  se_i <- .quadform_se(g_int, vcov_mat)
  ci_i <- l_int + c(-1, 1) * zcrit * se_i

  cbind(
    lab,
    data.frame(
      DDlog2MIC  = exp(l_int),
      SE_logDoD = se_i,
      CI_Lower  = exp(ci_i[1]),
      CI_Upper  = exp(ci_i[2]),
      P_value   = 2 * stats::pnorm(-abs(l_int / se_i)),
      stringsAsFactors = FALSE
    )
  )
}


# ---------- ADDITIVE (difference-of-differences on original scale) ----------
#' .difference_of_differences_additive (internal)
#' Additive Δ-of-Δs for 2x2 designs.
#' @keywords internal
#' @noRd
.difference_of_differences_additive <- function(mic, g_mic, vcov_mat, zcrit, mic_df) {
  lab <- .dod_label_from_mic_df(mic_df, mic_col = "MIC")
  if (is.null(lab) || length(mic) != 4L) return(.dod_empty_row("additive"))

  if (!is.matrix(g_mic) || nrow(g_mic) != 4L)
    stop("ordinalMIC: g_mic must be 4xp")
  if (!is.matrix(vcov_mat))
    stop("ordinalMIC: vcov_mat must be a square matrix")

  d_int <- (mic[4] - mic[2]) - (mic[3] - mic[1])
  g_int <- (g_mic[4,,drop=FALSE] - g_mic[2,,drop=FALSE]) -
    (g_mic[3,,drop=FALSE] - g_mic[1,,drop=FALSE])

  se_i <- .quadform_se(g_int, vcov_mat)
  ci_i <- d_int + c(-1, 1) * zcrit * se_i

  cbind(
    lab,
    data.frame(
      DDMIC = d_int,
      SE_DoD   = se_i,
      CI_Lower = ci_i[1],
      CI_Upper = ci_i[2],
      P_value  = 2 * stats::pnorm(-abs(d_int / se_i)),
      stringsAsFactors = FALSE
    )
  )
}




# DoD on ratio (log scale), across all 2×2 slices
.difference_of_differences_ratio_all <- function(log_mic, g_log_mic, vcov_mat, zcrit, mic_df) {
  slices <- .dod_enumerate_slices(mic_df, mic_col = "MIC")
  if (!length(slices)) return(.dod_empty_row("ratio"))

  rows <- lapply(slices, function(s) {
    i <- s$idx
    l_int <- (log_mic[i[4]] - log_mic[i[2]]) - (log_mic[i[3]] - log_mic[i[1]])
    g_int <- (g_log_mic[i[4], , drop=FALSE] - g_log_mic[i[2], , drop=FALSE]) -
      (g_log_mic[i[3], , drop=FALSE] - g_log_mic[i[1], , drop=FALSE])

    if (!all(is.finite(c(l_int, g_int)))) return(NULL)

    se_i <- .quadform_se(g_int, vcov_mat)
    ci_i <- l_int + c(-1,1)*zcrit*se_i

    data.frame(
      var1 = s$f1, var2 = s$f2,
      var1_lvlA = s$low1,  var1_lvlB = s$high1,
      var2_lvlC = s$low2,  var2_lvlD = s$high2,
      label = s$label,
      DDlog2MIC  = (l_int)/log(2),
      SE_logDoD = se_i/log(2),
      CI_Lower  = (ci_i[1])/log(2),
      CI_Upper  = (ci_i[2])/log(2),
      P_value   = 2*stats::pnorm(-abs(l_int/se_i)),
      stringsAsFactors = FALSE
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(.dod_empty_row("ratio"))
  do.call(rbind, rows)
}

# DoD on additive scale, across all 2×2 slices
.difference_of_differences_additive_all <- function(mic, g_mic, vcov_mat, zcrit, mic_df) {
  slices <- .dod_enumerate_slices(mic_df, mic_col = "MIC")
  if (!length(slices)) return(.dod_empty_row("additive"))

  rows <- lapply(slices, function(s) {
    i <- s$idx
    d_int <- (mic[i[4]] - mic[i[2]]) - (mic[i[3]] - mic[i[1]])
    g_int <- (g_mic[i[4], , drop=FALSE] - g_mic[i[2], , drop=FALSE]) -
      (g_mic[i[3], , drop=FALSE] - g_mic[i[1], , drop=FALSE])

    if (!all(is.finite(c(d_int, g_int)))) return(NULL)

    se_i <- .quadform_se(g_int, vcov_mat)
    ci_i <- d_int + c(-1,1)*zcrit*se_i

    data.frame(
      var1 = s$f1, var2 = s$f2,
      var1_lvlA = s$low1,  var1_lvlB = s$high1,
      var2_lvlC = s$low2,  var2_lvlD = s$high2,
      label = s$label,
      DDMIC = d_int,
      SE_DoD   = se_i,
      CI_Lower = ci_i[1],
      CI_Upper = ci_i[2],
      P_value  = 2*stats::pnorm(-abs(d_int/se_i)),
      stringsAsFactors = FALSE
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(.dod_empty_row("additive"))
  do.call(rbind, rows)
}

