# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

hardThresholdCpp <- function(u, thres) {
    .Call(`_TSLA_hardThresholdCpp`, u, thres)
}

shrinkGroupCpp <- function(u, g_idx) {
    .Call(`_TSLA_shrinkGroupCpp`, u, g_idx)
}

#' Calculate group norms
#'
#' Function to output group norms on the gamma coefficients
#' based on \code{g_idx} and \code{C2} matrix.
#'
#' @param u \code{C2*gamma.coef}, \code{gamma.coef} is the estimated node coefficient vector,
#' \code{C2} matrix is the output from
#' function \code{get_tree_object()}, which gives the weights of the groups.
#' @param g_idx Group structure matrix defined by the \code{C2} matrix.
#' See details in \code{get_tree_object()}.
#' @param type If \code{type == 1}, return sum of group norms; else return individual
#' norm for each group.
#' @return Sum of group norms or individual group norms.
#' @export
cal2norm <- function(u, g_idx, type) {
    .Call(`_TSLA_cal2norm`, u, g_idx, type)
}

SPGlinear <- function(y, X_1, X_2, C_1, C_1norm, C_2, C_2norm, g_idx, gamma_est, lambda, alpha, Mu, maxit, tol, verbose) {
    .Call(`_TSLA_SPGlinear`, y, X_1, X_2, C_1, C_1norm, C_2, C_2norm, g_idx, gamma_est, lambda, alpha, Mu, maxit, tol, verbose)
}

SPGlogistic <- function(y, X_1, X_2, C_1, C_1norm, C_2, C_2norm, g_idx, gamma_est, weight, lambda, alpha, Mu, maxit, tol, verbose) {
    .Call(`_TSLA_SPGlogistic`, y, X_1, X_2, C_1, C_1norm, C_2, C_2norm, g_idx, gamma_est, weight, lambda, alpha, Mu, maxit, tol, verbose)
}

warm_start_ls <- function(y, X_1, X_2, C_1, C_1norm, C_2, C_2norm, g_idx, gamma_init, lambda, alpha, Mu, maxit, tol, verbose) {
    .Call(`_TSLA_warm_start_ls`, y, X_1, X_2, C_1, C_1norm, C_2, C_2norm, g_idx, gamma_init, lambda, alpha, Mu, maxit, tol, verbose)
}

warm_start_logistic <- function(y, X_1, X_2, C_1, C_1norm, C_2, C_2norm, g_idx, gamma_init, weight, lambda, alpha, Mu, maxit, tol, verbose) {
    .Call(`_TSLA_warm_start_logistic`, y, X_1, X_2, C_1, C_1norm, C_2, C_2norm, g_idx, gamma_init, weight, lambda, alpha, Mu, maxit, tol, verbose)
}

