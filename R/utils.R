#' #' Make a concise DoD label from two group level vectors
#' #' @param g1,g2 named character vectors of factor levels
#' #' @return single character string, eg  "strain: Mut-WT | treatment: Salt–Mock"
#' #' @noRd
#' .make_dod_label <- function(g1, g2) {
#'   vars <- names(g1)
#'   dif  <- vars[g1 != g2]
#'   if (length(dif) == 0) {
#'     "no difference"
#'   } else {
#'     paste(
#'       sprintf("%s: %s\u2013%s", dif, g2[dif], g1[dif]),   # em-dash
#'       collapse = " | ")
#'   }
#' }

#' Build a readable label for a 2×2 difference-of-differences
#'
#' @param groups character vector length 4, in the order:
#'   1 = (low1, low2), 2 = (high1, low2),
#'   3 = (low1, high2), 4 = (high1, high2)
#' @param sep    separator used inside each group string (default ":")
#' @return a single character string
#' @noRd
# .dod_make_label <- function(groups, sep = ":") {
#   if (length(groups) != 4L) return(NA_character_)
#
#   # split "WT:None" → c("WT","None"), keep the original variable names if present
#   split_g <- strsplit(groups, sep, fixed = TRUE)
#   nvars   <- length(split_g[[1L]])
#   if (!all(vapply(split_g, length, integer(1)) == nvars))
#     stop("All group strings must have the same number of factors.")
#
#   # if names are present in mic_solve$newdata, use them; else V1, V2, …
#   vars <- names(split_g[[1L]])
#   if (is.null(vars)) vars <- paste0("V", seq_len(nvars))
#
#   g1 <- split_g[[1L]]     # baseline corner
#   g4 <- split_g[[4L]]     # opposite corner
#   names(g1) <- names(g4) <- vars
#
#   diff_vars <- vars[g1 != g4]
#
#   if (length(diff_vars) == 0) {
#     "no difference"
#   } else if (length(diff_vars) == 1) {
#     v <- diff_vars
#     sprintf("%s: %s \u2013 %s", v, g4[[v]], g1[[v]])   # en-dash
#   } else {
#     # two or more factors changed → show full group strings
#     sprintf("%s \u2013 %s", groups[1L], groups[4L])
#   }
# }


#' @noRd
.dod_all_base <- function(mic_tbl, mic_col = "MIC") {

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
            label = sprintf("%s-%s vs. %s-%s",
                            #v1,
                            lv1[a],
                            lv1[b],
                            #v2,
                            lv2[c],
                            lv2[d]),
            stringsAsFactors = FALSE
          )
        }
    }
  }

  if (k == 0L) stop("No 2x2 sub-grids found that have >=2 levels per factor.")
  do.call(rbind, out)
}
## ---------------------------------------------------------------
# is_design <- vapply(object$mic_estimates,     # logical vector
#                     function(col) !(is.numeric(col) || is.integer(col)),
#                     logical(1))

#vars     <- names(object$mic_estimates)[is_design]      # e.g. "strain","treatment"

#
#
# tbl  <- object$dod_ratio_results          # usually one row
# vars <- names(object$newdata)             # e.g. c("strain","treatment")
#
# mic_tbl  <- object$mic_estimates[c(vars, "MIC")]
#
# dod_tbl <- dod_all_base(mic_tbl, mic_col = "MIC")
#
# tbl <- merge(tbl, dod_tbl, by = c("f1", "f2", "low1", "high1", "low2", "high2"),
#              all.x = TRUE, sort = FALSE)
