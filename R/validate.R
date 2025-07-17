#' .validate_inputs (internal)
#'
#' Ensures the model, `newdata`, and concentration variable are supplied in a
#' minimally sane form before the heavy lifting begins.  The concentration
#' column does **not** have to appear in `newdata`; it will be injected later
#' with a placeholder.  If the column *is* present, we simply check that it's
#' numeric.
#'
#' @keywords internal
#' @noRd
.validate_inputs <- function(clm_fit, newdata, conc_name) {
  if (!requireNamespace("ordinal", quietly = TRUE)) {
    stop("Package 'ordinal' is required but not installed.")
  }

  if (!inherits(clm_fit, "clm")) {
    stop("'clm_fit' must be created with ordinal::clm().")
  }

  if (!is.data.frame(newdata)) {
    stop("'newdata' must be a data frame.")
  }

  # Allow the concentration column to be absent here; placeholder is added
  # later.  If the user *does* include it, enforce that it's numeric.
  if (conc_name %in% names(newdata)) {
    if (!is.numeric(newdata[[conc_name]])) {
      stop(sprintf("Column '%s' must be numeric if supplied in 'newdata'.",
                   conc_name))
    }
  }
  invisible(TRUE)
}
