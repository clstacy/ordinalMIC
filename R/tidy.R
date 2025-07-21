#' tidy method for mic_solve objects
#' @importFrom tibble tibble
#' @importFrom tibble as_tibble
#' @param x An object of class `mic_solve`.
#' @param table A character string specifying which table to return.
#'   Options are `"mic"`, `"delta"`, `"ratio"`, `"DoD_delta"`, and `"DoD_ratio"`.
#' @param ... Additional arguments (not used).
#' @return A tibble containing the requested table from the `mic_solve` object.
#' @export
tidy.mic_solve <- function(x, table = c("mic", "delta", "ratio",
                                        "DoD_delta", "DoD_ratio"), ...) {
  table <- match.arg(table)
  out <- switch(table,
                mic       = x$mic_estimates,
                delta     = x$delta_mic_results,
                ratio     = x$ratio_mic_results,
                DoD_delta = x$dod_delta_results,
                DoD_ratio = x$dod_ratio_results)
  tibble::as_tibble(out)
}
