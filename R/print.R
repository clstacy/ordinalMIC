#' @export
print.mic_solve <- function(x, ...) {
  cat(cli::rule(left = "ordinalMIC summary"), "\n")
  cat("Call  : ", deparse(x$call), "\n\n")

  g <- nrow(x$mic_estimates)
  cat("Groups:", g,   "(", paste(names(x$newdata), collapse = " x "), ")\n")
  cat("Alpha : ", attr(x$model$alpha, "alpha"), "\n")
  cat("Link  : ", if (identical(x$transform_fun, identity)) "identity"
      else deparse(x$transform_fun), "\n\n")

  cli::cat_line(cli::style_bold("MIC table"))
  print(utils::head(x$mic_estimates, 6), row.names = FALSE)
  cat("\n")
  cli::cat_line(cli::style_bold("deltaMIC table"))
  print(utils::head(x$delta_mic_results, 6), row.names = FALSE)
  cat("\n")
  cli::cat_line(cli::style_bold("FC MIC table"))
  print(utils::head(x$ratio_mic_results, 6), row.names = FALSE)
  if (!is.null(x$dod_delta_results)) {
    cli::cat_line(cli::style_bold("DoD deltaMIC table"))
    print(utils::head(x$dod_delta_results, 6), row.names = FALSE)
    cat("\n")
  }
  if (!is.null(x$dod_ratio_results)) {
    cli::cat_line(cli::style_bold("DoD ratioMIC table"))
    print(utils::head(x$dod_ratio_results, 6), row.names = FALSE)
    cat("\n")
  }

  invisible(x)
}
