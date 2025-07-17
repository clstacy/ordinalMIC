#' @importFrom tibble tibble
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
