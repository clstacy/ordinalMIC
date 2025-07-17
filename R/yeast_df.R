#' Simulated Yeast Stress-Survival Data
#'
#' An example dataset with ordinal spot scores (0 = no growth … 4 = full
#' growth) for treated and untreated wild-type and mutant yeast
#' across a peroxide gradient.
#'
#' @format A data frame with 132 rows and 4 variables:
#' \describe{
#'   \item{score}{Ordered factor \code{0 < 1 < 2 < 3 < 4}}
#'   \item{conc}{Numeric concentration units}
#'   \item{strain}{Factor: \code{"WT"} or \code{"Mut"}}
#'   \item{treatment}{Factor: \code{"None"} or \code{"Salt"}}
#' }
#'
#' @usage data(yeast_df)
#'
#' @examples
#' data(yeast_df)
#' head(yeast_df)
"yeast_df"
