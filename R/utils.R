#' .dod_all_base (internal)
#' Ratio-of-ratios (multiplicative DoD) for 2×2 designs.
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
.dod_all_base <- function(mic_tbl, mic_col = "MIC") {

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



