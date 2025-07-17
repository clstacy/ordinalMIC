#' .build_model_matrix (internal)
#' @keywords internal
#' @noRd
.build_model_matrix <- function(clm_fit, newdata, conc_name,
                                placeholder_val) {
  mm_terms <- stats::delete.response(stats::terms(clm_fit))
  tmp      <- newdata
  tmp[[conc_name]] <- placeholder_val
  stats::model.matrix(mm_terms, data = tmp, xlev = clm_fit$xlevels)
}
