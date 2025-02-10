#' Synthesic for the regression example
#'
#' Synthetic data used to illustrate how to use TSLA with regression.
#'
#' @name RegressionExample
#' @docType data
#' @usage data(RegressionExample)
#' @keywords data
#' @format List containing the following elements:
#' \describe{
#'   \item{tree.org}{Original tree structure with 42 leaf nodes and 5 different levels.}
#'   \item{x.org}{Original design matrix with 42 binary features and 400 observations.}
#'   \item{y}{Continuous response of length 400.}
#' }
"RegressionExample"


#' Synthesic for the classification example
#'
#' Synthetic data used to illustrate how to use TSLA with classification.
#'
#' @name ClassificationExample
#' @docType data
#' @usage data(ClassificationExample)
#' @keywords data
#' @format List containing the following elements:
#' \describe{
#'   \item{tree.org}{Original tree structure with 42 leaf nodes and 5 different levels.}
#'   \item{x.org}{Original design matrix with 42 binary features and 400 observations.}
#'   \item{x1}{Unpenalized covariate.}
#'   \item{y}{Continuous response of length 400.}
#' }
"ClassificationExample"
