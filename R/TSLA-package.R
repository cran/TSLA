#' Tree-Guided Rare Feature Selection and Logic Aggregation
#'
#' This package provides functions and visualization tools for fitting the
#' Tree-Guided Rare Feature Selection and Logic Aggregation model.
#'
#' @author Jianmin Chen \email{jianminc000@gmail.com}, Kun Chen
#' @name TSLA-package
#' @docType package

#' @useDynLib TSLA, .registration=TRUE
#' @import Rcpp
#' @import RcppArmadillo
#' @import Matrix
#' @importFrom pROC auc roc coords
#' @importFrom PRROC pr.curve
#' @importFrom stats quantile
#' @importFrom stats sd
#' @importFrom stats rnorm
#' @importFrom stats rbinom
#' @importFrom stats var
#' @importFrom stats as.formula
#' @importFrom stats model.matrix
#' @importFrom ape as.phylo
#' @importFrom ape tiplabels
#' @importFrom ape nodelabels
#' @importFrom phytools phylo.to.map
#' @importFrom data.tree as.Node
NULL
