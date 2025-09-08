#' .dod_all_multiplicative (internal)
#' Ratio-of-ratios (multiplicative DoD) for 2x2 designs.
#' @param mic_tbl  Data frame with MIC estimates and design factors.
#' @param mic_col  Name of the MIC column in `mic_tbl` (default `"MIC"`).
#' @return A data frame with columns:
#'   * `f1`, `f2` – factor names of the two design factors
#'   * `low1`, `high1` – levels of the first factor
#'   * `low2`, `high2` – levels of the second factor
#'   * `ratio` – ratio of ratios (DoD)
#'   * `log2_ratio` – log2 of the ratio
#'   * `label` – a human-readable label for the comparison
#' @keywords internal
#' @noRd
.dod_all_multiplicative <- function(mic_tbl, mic_col = "MIC") {

  # remove SE and CI columns from mic table
  mic_tbl <- mic_tbl[, c(mic_col, setdiff(names(mic_tbl), c(mic_col, "SE_LP", "CI_Lower", "CI_Upper")))]
  vars <- setdiff(names(mic_tbl), mic_col)   # every design column
  if (length(vars) < 2L)
    stop("Need at least two design factors for a DoD.")

  mic_vals <- mic_tbl[[mic_col]]

  # helper: fetch MIC for a combination held in named character vector
  get_mic <- function(levels) {
    idx <- rep(TRUE, nrow(mic_tbl))
    for (v in names(levels))
      idx <- idx & (as.character(mic_tbl[[v]]) == levels[[v]])
    mic_vals[idx][1L]          # assume unique match
  }

  out <- list()                # will store rows, then rbind
  k   <- 0L

  for (i in seq_along(vars)[-length(vars)]) {
    for (j in (i + 1L):length(vars)) {

      v1 <- vars[i];  lv1 <- unique(as.character(mic_tbl[[v1]]))
      v2 <- vars[j];  lv2 <- unique(as.character(mic_tbl[[v2]]))

      if (length(lv1) < 2L || length(lv2) < 2L) next    # need 2 levels each

      # other factors: hold at *first* observed level
      baseline <- lapply(setNames(vars, vars),
                         function(v) as.character(mic_tbl[[v]][1L]))

      # iterate over every pair of levels within v1 and v2
      for (a in 1:(length(lv1) - 1L)) for (b in (a + 1L):length(lv1))
        for (c in 1:(length(lv2) - 1L)) for (d in (c + 1L):length(lv2)) {

          comb <- baseline

          comb[[v1]] <- lv1[a]; comb[[v2]] <- lv2[c]   # low1 / low2
          a_i <- get_mic(comb)

          comb[[v1]] <- lv1[b]                         # high1  / low2
          b_i <- get_mic(comb)

          comb[[v1]] <- lv1[a]; comb[[v2]] <- lv2[d]   # low1  / high2
          a_j <- get_mic(comb)

          comb[[v1]] <- lv1[b];                        # high1 / high2
          b_j <- get_mic(comb)

          ratio <- (b_j / a_j) / (b_i / a_i)

          k <- k + 1L
          out[[k]] <- data.frame(
            f1    = v1,  f2 = v2,
            low1  = lv1[a], high1 = lv1[b],
            low2  = lv2[c], high2 = lv2[d],
            ratio = ratio,
            log2_ratio = log2(ratio),
            label = sprintf("%s - %s vs. %s - %s",
                            #v1,
                            lv1[b],
                            lv1[a],
                            #v2,
                            lv2[d],
                            lv2[c]),
            stringsAsFactors = FALSE
          )
        }
    }
  }

  if (k == 0L) stop("No 2x2 sub-grids found that have >=2 levels per factor.")
  do.call(rbind, out)
}


#' .dod_all_additive (internal)
#' Additive Δ-of-Δs for 2x2 designs.
#' @param mic_tbl  Data frame with MIC estimates and design factors.
#' @param mic_col  Name of the MIC column in `mic_tbl` (default `"MIC"`).
#' @return A data frame with columns:
#' * `f1`, `f2` – factor names of the two design factors
#' * `low1`, `high1` – levels of the first factor
#' * `low2`, `high2` – levels of the second factor
#' * `delta` – difference of differences (DoD)
#' * `label` – a human-readable label for the comparison
#' @keywords internal
#' @noRd
.dod_all_additive <- function(mic_tbl, mic_col = "MIC") {

  # remove SE and CI columns from mic table
  mic_tbl <- mic_tbl[, c(mic_col, setdiff(names(mic_tbl), c(mic_col, "SE_LP", "CI_Lower", "CI_Upper")))]
  vars <- setdiff(names(mic_tbl), mic_col)   # every design column
  if (length(vars) < 2L)
    stop("Need at least two design factors for a DoD.")

  mic_vals <- mic_tbl[[mic_col]]

  # helper: fetch MIC for a combination held in named character vector
  get_mic <- function(levels) {
    idx <- rep(TRUE, nrow(mic_tbl))
    for (v in names(levels))
      idx <- idx & (as.character(mic_tbl[[v]]) == levels[[v]])
    mic_vals[idx][1L]          # assume unique match
  }

  out <- list()                # will store rows, then rbind
  k   <- 0L

  for (i in seq_along(vars)[-length(vars)]) {
    for (j in (i + 1L):length(vars)) {

      v1 <- vars[i];  lv1 <- unique(as.character(mic_tbl[[v1]]))
      v2 <- vars[j];  lv2 <- unique(as.character(mic_tbl[[v2]]))

      if (length(lv1) < 2L || length(lv2) < 2L) next    # need 2 levels each

      # other factors: hold at *first* observed level
      baseline <- lapply(setNames(vars, vars),
                         function(v) as.character(mic_tbl[[v]][1L]))

      # iterate over every pair of levels within v1 and v2
      for (a in 1:(length(lv1) - 1L)) for (b in (a + 1L):length(lv1))
        for (c in 1:(length(lv2) - 1L)) for (d in (c + 1L):length(lv2)) {

          comb <- baseline

          comb[[v1]] <- lv1[a]; comb[[v2]] <- lv2[c]   # low1 / low2
          a_i <- get_mic(comb)


          comb[[v1]] <- lv1[b]                         # high1  / low2
          b_i <- get_mic(comb)

          comb[[v1]] <- lv1[a]; comb[[v2]] <- lv2[d]   # low1  / high2
          a_j <- get_mic(comb)

          comb[[v1]] <- lv1[b];                        # high1 / high2
          b_j <- get_mic(comb)

          delta <- (b_j - a_j) - (b_i - a_i)

          k <- k + 1L
          out[[k]] <- data.frame(
            f1    = v1,  f2 = v2,
            low1  = lv1[a], high1 = lv1[b],
            low2  = lv2[c], high2 = lv2[d],
            delta = delta,
            label = sprintf("%s - %s vs. %s - %s",
                            #v1,
                            lv1[b],
                            lv1[a],
                            #v2,
                            lv2[d],
                            lv2[c]),
            stringsAsFactors = FALSE
          )
        }
    }
  }

  if (k == 0L) stop("No 2x2 sub-grids found that have >=2 levels per factor.")
  do.call(rbind, out)
}

#
#
# # Map a 1xp gradient row to the full coefficient space of V (vcov)
# # by name, zero-filling any missing coefficients, and reordering to V's cols.
# .quadform_se <- function(g_row, V) {
#   if (is.null(colnames(V))) stop("ordinalMIC: vcov_mat must have colnames")
#   # coerce g_row to a named numeric vector
#   if (is.matrix(g_row)) {
#     gi <- drop(g_row)
#     names(gi) <- colnames(g_row)
#   } else {
#     gi <- as.numeric(g_row)
#     if (is.null(names(gi))) stop("ordinalMIC: gradient row must have names")
#   }
#   full <- numeric(ncol(V))
#   names(full) <- colnames(V)
#   common <- intersect(names(gi), names(full))
#   if (length(common) == 0L) stop("ordinalMIC: no overlapping parameter names between gradient and vcov")
#   full[common] <- gi[common]
#   sqrt(drop(t(full) %*% V %*% full))
# }
#
#
# # Align a (row) gradient to the column space of V (vcov), by NAMES.
# # - Accepts numeric vector or 1xp matrix; returns a 1xncol(V) matrix.
# .align_grad_row_to_vcov <- function(g_row, V) {
#   if (!is.matrix(V)) stop("ordinalMIC: vcov_mat must be a matrix")
#   vcols <- colnames(V)
#   if (is.null(vcols)) stop("ordinalMIC: vcov_mat must have colnames")
#
#   # coerce g_row -> named numeric
#   if (is.matrix(g_row)) {
#     gi <- drop(g_row)               # length p
#     names(gi) <- colnames(g_row)
#   } else {
#     gi <- as.numeric(g_row)
#     if (is.null(names(gi))) stop("ordinalMIC: gradient row must be named")
#   }
#
#   out <- matrix(0, nrow = 1, ncol = ncol(V))
#   colnames(out) <- vcols
#
#   if (length(gi)) {
#     common <- intersect(names(gi), vcols)
#     if (length(common) == 0L) {
#       stop("ordinalMIC: no overlapping parameter names between gradient and vcov")
#     }
#     out[, common] <- gi[common]
#   }
#   out
# }
#
# # Align a stack of gradients (mxp) to V, by names; returns mxncol(V).
# .align_grad_mat_to_vcov <- function(G, V) {
#   if (is.null(colnames(G))) {
#     stop("ordinalMIC: gradient matrix must have colnames matching model coefficients")
#   }
#   vcols <- colnames(V)
#   if (is.null(vcols)) stop("ordinalMIC: vcov_mat must have colnames")
#
#   out <- matrix(0, nrow = nrow(G), ncol = ncol(V))
#   colnames(out) <- vcols
#   rownames(out) <- rownames(G)
#
#   common <- intersect(colnames(G), vcols)
#   if (length(common) == 0L) stop("ordinalMIC: gradient and vcov share no coefficient names")
#
#   out[, common, drop = FALSE] <- G[, common, drop = FALSE]
#   out
# }
#
# # Quadratic-form SE for a single row gradient using name-alignment.
# .quadform_se <- function(g_row, V) {
#   gA <- .align_grad_row_to_vcov(g_row, V)  # 1xk
#   sqrt(drop(gA %*% V %*% t(gA)))
# }

# Align a 1xp gradient row (or a named vector) to V by names; return 1xk matrix
.align_grad_row_to_vcov <- function(g_row, V) {
  if (!is.matrix(V)) stop("ordinalMIC: vcov_mat must be a matrix")
  vcols <- colnames(V)
  if (is.null(vcols)) stop("ordinalMIC: vcov_mat must have colnames")

  # Coerce to a 1xp matrix with colnames
  if (!is.matrix(g_row)) {
    if (is.null(names(g_row))) stop("ordinalMIC: gradient row must be named")
    g_row <- matrix(as.numeric(g_row), nrow = 1L)
    colnames(g_row) <- names(g_row)
  }

  Gcols <- colnames(g_row)
  if (is.null(Gcols)) stop("ordinalMIC: gradient row must have colnames")

  out <- matrix(0, nrow = 1L, ncol = ncol(V))
  colnames(out) <- vcols

  common <- intersect(Gcols, vcols)
  if (length(common) == 0L) stop("ordinalMIC: no overlapping parameter names between gradient and vcov")

  vi <- match(common, vcols)
  gi <- match(common, Gcols)
  out[, vi] <- g_row[, gi, drop = FALSE]
  out
}

# Align an mxp gradient matrix (or a named vector) to V by names; return mxk
.align_grad_mat_to_vcov <- function(G, V) {
  if (!is.matrix(V)) stop("ordinalMIC: vcov_mat must be a matrix")
  vcols <- colnames(V)
  if (is.null(vcols)) stop("ordinalMIC: vcov_mat must have colnames")

  # Coerce vectors to 1xp matrices
  if (!is.matrix(G)) {
    if (is.null(names(G))) stop("ordinalMIC: gradient matrix must have colnames")
    G <- matrix(as.numeric(G), nrow = 1L)
    colnames(G) <- names(G)
  }

  Gcols <- colnames(G)
  if (is.null(Gcols)) stop("ordinalMIC: gradient matrix must have colnames")

  out <- matrix(0, nrow = nrow(G), ncol = ncol(V))
  rownames(out) <- rownames(G)
  colnames(out) <- vcols

  common <- intersect(Gcols, vcols)
  if (length(common) == 0L) stop("ordinalMIC: gradient and vcov share no coefficient names")

  vi <- match(common, vcols)
  gi <- match(common, Gcols)

  # single block assignment; both are guaranteed 2-D now
  out[, vi] <- G[, gi, drop = FALSE]
  out
}

# Quadratic-form SE for a single (possibly vector) gradient row
.quadform_se <- function(g_row, V) {
  gA <- .align_grad_row_to_vcov(g_row, V)  # 1xk
  sqrt(drop(gA %*% V %*% t(gA)))
}


# Enumerate all 2x2 sub-grids (indices) across factor pairs & level pairs.
# Returns a list of records with $idx (length 4), factor names/levels, and a label.
.dod_enumerate_slices <- function(mic_df, mic_col = "MIC") {
  design_cols <- setdiff(names(mic_df), c(mic_col, "SE_LP","CI_Lower","CI_Upper"))
  if (length(design_cols) < 2L) return(list())

  as_chr <- function(x) if (is.factor(x)) as.character(x) else as.character(x)
  DF <- mic_df
  for (v in design_cols) DF[[v]] <- as_chr(DF[[v]])

  out <- list(); k <- 0L

  for (i in seq_len(length(design_cols) - 1L)) {
    for (j in (i + 1L):length(design_cols)) {
      f1 <- design_cols[i]; lv1 <- unique(DF[[f1]])
      f2 <- design_cols[j]; lv2 <- unique(DF[[f2]])
      if (length(lv1) < 2L || length(lv2) < 2L) next

      # Hold all *other* factors fixed at their first observed level
      others <- setdiff(design_cols, c(f1, f2))
      baseline <- setNames(lapply(others, function(v) DF[[v]][1L]), others)

      # iterate over all pairs of levels within f1 and f2
      for (a in 1:(length(lv1)-1L)) for (b in (a+1L):length(lv1))
        for (c in 1:(length(lv2)-1L)) for (d in (c+1L):length(lv2)) {

          sel <- function(v,val) DF[[v]] == val
          mask_base <- if (length(others)) {
            Reduce("&", lapply(others, function(v) DF[[v]] == baseline[[v]]))
          } else rep(TRUE, nrow(DF))

          idx_ac <- which(sel(f1, lv1[a]) & sel(f2, lv2[c]) & mask_base)
          idx_ad <- which(sel(f1, lv1[a]) & sel(f2, lv2[d]) & mask_base)
          idx_bc <- which(sel(f1, lv1[b]) & sel(f2, lv2[c]) & mask_base)
          idx_bd <- which(sel(f1, lv1[b]) & sel(f2, lv2[d]) & mask_base)

          if (length(idx_ac)==1L && length(idx_ad)==1L &&
              length(idx_bc)==1L && length(idx_bd)==1L) {
            k <- k + 1L
            out[[k]] <- list(
              f1=f1, f2=f2,
              low1=lv1[a], high1=lv1[b],
              low2=lv2[c], high2=lv2[d],
              idx=c(idx_ac, idx_ad, idx_bc, idx_bd),  # canonical order: (a,c),(a,d),(b,c),(b,d)
              label=sprintf("%s: %s vs %s  x  %s: %s vs %s",
                            f1, lv1[b], lv1[a], f2, lv2[d], lv2[c])
            )
          }
        }
    }
  }
  out
}


if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "._bar_pos_", "y_obs", "y_cens", "x_pos", "y", "color_val", "cens", "facet",
    "log2Ratio_MIC", "DDMIC"
  ))
}

#' Internal: relevel a factor only if the reference level exists
#'
#' @param x   A factor or character vector.
#' @param ref A single character string naming the desired reference level.
#' @return A factor with the same levels (re-leveled if `ref` exists).
#' @keywords internal
#' @noRd
.relevel_if_present <- function(x, ref) {
  x <- as.factor(x)
  if (ref %in% levels(x)) stats::relevel(x, ref = ref) else x
}

#' Internal imports
#' @keywords internal
#' @importFrom stats na.omit
#' @importFrom utils globalVariables
#' @noRd
NULL

