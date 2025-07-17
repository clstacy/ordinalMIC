#' @export
print.mic_solve <- function(x, ...) {
  cat(cli::rule(left = "ordinalMIC summary"), "\n")
  cat("Call  : ", deparse(x$call), "\n\n")

  g <- nrow(x$mic_estimates)
  cat("Groups:", g,   "(", paste(names(x$newdata), collapse = " × "), ")\n")
  cat("Alpha : ", attr(x$mic_estimates, "alpha"), "\n")
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

  invisible(x)
}
