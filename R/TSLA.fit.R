#Rcpp::sourceCpp('./SPGlinear.cpp')
#Rcpp::sourceCpp('./SPGlogistic.cpp')
#source('./get_tree_object.R')

### internal function
###===========================================================================
#### control the model fitting: set convergence parameter
TSLA.control <- function(tol = 1e-5, maxit = 10000L, mu = 1e-3, verbose = FALSE) {
  list(tol = tol, maxit = maxit, mu = mu, verbose = verbose)
}

#### control the model fitting: set tuning parameter options
TSLA.modstr <- function(lambda = NULL, lambda.min.ratio = 1e-4,
                        nlambda = 50, alpha = seq(0, 1, length.out = 10)) {
  list(lambda = lambda,
       lambda.min.ratio = lambda.min.ratio,
       nlambda = nlambda,
       alpha = sort(alpha)) # make sure alpha values are in ascending order
}

### internal function
###===========================================================================
### function to estimate an upper bound for lambda
estlambda <- function(y, X, A, g_idx, wg, wj){
  # X: expanded design matrix for binary rare features
  n.g <- nrow(g_idx)
  X.gamma <- X %*% A
  gamma.g <- rep(0, n.g)
  for(i in 1:n.g){
    X.g <- X.gamma[, g_idx[i, 1]:g_idx[i, 2]]
    gamma.g[i] <- sqrt(sum((t(X.g) %*% y)^2))/wg[i]
  }
  gamma.g <- max(gamma.g)
  p <- ncol(X)
  gamma.j <- rep(0, p)
  for(j in 1:p){
    xj <- X[, j]
    gamma.j[j] <- abs(t(xj) %*% y)/wj[j]
  }
  gamma.j <- max(gamma.j)
  return(max(gamma.g, gamma.j))
}

### internal function
###===========================================================================
## Solve the optimization problem in TSLA with a specific alpha and a sequence of lambda
## Used to call the Rcpp function with warm start
## Return estimation of the gamma coefficients
TSLA <- function(y, X_1 = NULL, X_2, C_1, C_1norm, C_2, C_2norm,
                 g_idx, A, family = c('ls', 'logit'),
                 gamma.init = NULL, weight = NULL, lambda = 1, alpha = 0,
                 control = list()){
  # get control parameters
  maxit <- control$maxit
  mu <- control$mu
  tol <- control$tol
  verbose <- control$verbose

  # add intercept column in X1
  if(is.null(X_1)){
    X_1 <- matrix(1, nrow = nrow(as.matrix(y)), ncol = 1)
  }else{
    X_1 <- cbind(rep(1, length(y)), X_1)
  }

  # set default for initial value: 0 vector
  if(is.null(gamma.init)){
    gamma.init <- rep(0, ncol(X_1) + ncol(A))
  }

  # call Rcpp function with warm start based on family
  family <- match.arg(family)
  fit <- switch(family,
                ls = warm_start_ls(y, X_1, X_2 %*% A, C_1, C_1norm, C_2, C_2norm,
                                     g_idx, gamma.init, lambda, alpha, mu, maxit, tol, verbose),
                logit = warm_start_logistic(y, X_1, X_2 %*% A, C_1, C_1norm, C_2, C_2norm,
                                    g_idx, gamma.init, weight, lambda, alpha, mu, maxit, tol, verbose))
                #ls = SPGlinear(y, X_1, X_2 %*% A, C_1, C_1norm, C_2, C_2norm,
                #     g_idx, gamma.init, lambda, alpha, mu, maxit, tol, verbose),
                #logit = SPGlogistic(y, X_1, X_2 %*% A, C_1, C_1norm, C_2, C_2norm,
                #              g_idx, gamma.init, weight, lambda, alpha, mu, maxit, tol, verbose))
  return(list(gamma.est = fit))

}


### Fit the regularization paths for TSLA learning problem for a sequence of alpha and lambda values
### Containing steps from reparameterization, and regularization
#' Solve the TSLA optimization problem
#'
#' Find the solutions with a Smoothing Proximal Gradient (SPG) algorithm
#' for a sequence of \eqn{\alpha} and \eqn{\lambda} values.
#'
#' We adopt the warm start technique to speed up the calculation.
#' The warm start is applied with a fixed value of \eqn{\alpha} and a
#' descending sequence of \eqn{\lambda}.
#'
#' The objective function for "ls" is
#' \deqn{1/2 RSS+\lambda(\alpha P(\beta)+(1-\alpha) P(\gamma)),}
#' subject to \eqn{\beta=A\gamma}.
#' The objective function for "logit" is
#' \deqn{-loglik + \lambda(\alpha P(\beta)+(1-\alpha) P(\gamma)),}
#' subject to \eqn{\beta=A\gamma}. Note that, in this package, the input parameter "alpha" is the
#' tuning parameter for the generalized lasso penalty.
#'
#' Details for "penalty" option:
#'
#' For \code{penalty = "CL2"}, see details for the
#' "Child-l2" penalty in the main paper.
#'
#' For \code{penalty = "RFS-Sum"}, the theoretical optimal weights are used.
#' Please check the details in paper
#' "Rare feature selection in high dimensions".
#'
#' @param y Response in matrix form, continuous for \code{family = "ls"} and binary (0/1)
#' for \code{family = "logit"}.
#' @param X_1 Design matrix for unpenalized features (excluding intercept). Need to be in the matrix form.
#' @param X_2 Expanded design matrix for \code{penalty = "CL2"}; Original design matrix
#' for \code{penalty = "RFS-Sum"}. Need to be in the matrix form.
#' @param treemat Expanded tree structure in matrix form for
#' \code{penalty = "CL2"}; Original structure for \code{penalty = "RFS-Sum"}.
#' @param family Two options. Use "ls" for least square problems and "logit"
#' for logistic regression problems.
#' @param penalty Two options for group penalty on \eqn{\gamma}, "CL2" or "RFS-Sum".
#' @param gamma.init Initial value for the optimization. Default is a zero vector.
#' The length should equal to 1+\code{ncol(X_1)}+\code{ncol(A)}.
#' See details of A in \code{get_tree_obj()}.
#' @param weight A vector of length two and it is used for logistic regression
#' only. The first element corresponds to weight of y=1 and the
#' second element corresponds to weight of y=0.
#' @param group.weight User-defined weights for group penalty. Need to be a vector
#' and the length equals to the number of groups.
#' @param feature.weight User-defined weights for each predictor after expansion.
#' @param control A list of parameters controlling algorithm convergence. Default values:
#' \code{tol = 1e-5}, convergence tolerance; \code{maxit = 10000},
#' maximum number of iterations; \code{mu = 1e-3}, smoothness parameter in SPG.
#' @param modstr A list of parameters controlling tuning parameters. Default values:
#' \code{lambda = NULL}. If lambda is not provided, the package will give a default lambda sequence;
#' \code{lambda.min.ratio = 1e-04}, smallest value for lambda as
#' a fraction of lambda.max (given by default when lambda is NULL);
#' \code{nlambda = 50},
#' number of lambda values (equal spacing on log scale) used when lambda is NULL;
#' \code{alpha = seq(0, 1, length.out = 10)}, sequence of alpha. Here, alpha is
#' tuning parameter for generalized lasso penalty and 1-alpha is the tuning
#' parameter for group lasso penalty.
#'
#' @return A list of model fitting results.
#' \item{gammacoef}{Estimation for \eqn{\gamma}.}
#' \item{groupnorm}{Weighted norms for each group.}
#' \item{lambda.seq}{Sequence of \eqn{\lambda} values.}
#' \item{alpha.seq}{Tuning parameter sequence for the generalized lasso penalty.}
#' \item{rmid}{Column index for all zero features.}
#' \item{family}{Option of \code{family}.}
#' \item{cov.name}{Names for unpenalized features.}
#' \item{bin.name}{Names for binary feautres.}
#' \item{tree.object}{Outputs from \code{get_tree_obj()}.}
#'
#' @references
#' Chen, J., Aseltine, R. H., Wang, F., & Chen, K. (2024).
#' \emph{Tree-Guided Rare Feature Selection and Logic Aggregation with
#' Electronic Health Records Data. Journal of the American Statistical Association 119(547), 1765-1777},
#' \doi{10.1080/01621459.2024.2326621}.\cr
#' Chen, X., Q. Lin, S. Kim, J. G. Carbonell, and E. P. Xing (2012).
#' \emph{Smoothing proximal gradient method for general structured sparse regression.
#' The Annals of Applied Statistics 6(2), 719–752},
#' \doi{10.1214/11-AOAS514}.\cr
#' Yan, X. and J. Bien (2021).
#' \emph{Rare feature selection in high dimensions.
#' Journal of the American Statistical Association 116(534), 887–900},
#' \doi{10.1080/01621459.2020.1796677}.\cr
#' @export
#'
#' @examples
#' # Load the synthetic data
#' data(RegressionExample)
#'
#' tree.org <- RegressionExample$tree.org   # original tree structure
#' x2.org <- RegressionExample$x.org      # original design matrix
#' y <- RegressionExample$y            # response
#'
#' # Do the tree-guided expansion
#' expand.data <- getetmat(tree.org, x2.org)
#' x2 <- expand.data$x.expand              # expanded design matrix
#' tree.expand <- expand.data$tree.expand  # expanded tree structure
#'
#' # specify some model parameters
#' set.seed(100)
#' control <- list(maxit = 100, mu = 1e-3, tol = 1e-5, verbose = FALSE)
#' # fit model with a pair of lambda and alpha
#' modstr <- list(lambda = 1,  alpha = 0.1)
#' x1 <- NULL
#' fit1 <- TSLA.fit(y, x1, x2, tree.expand, family = 'ls',
#'                  penalty = 'CL2',
#'                  gamma.init = NULL, weight = NULL,
#'                  group.weight = NULL, feature.weight = NULL,
#'                  control, modstr)
#' # get group norms from fit1
#' fit1$groupnorm



TSLA.fit <- function(y, X_1 = NULL, X_2, treemat, family = c('ls', 'logit'),
                     penalty = c('CL2', 'RFS-Sum'),
                     gamma.init = NULL, weight = NULL,
                     group.weight = NULL,
                     feature.weight = NULL,
                     control = list(), modstr = list()){
  ## X_1: design matrix for unpenalized covariates
  ## X_2: expanded design matrix for binary rare features

  control <- do.call("TSLA.control", control)
  modstr <- do.call("TSLA.modstr", modstr)

  # record function call
  Call <- match.call()
  family <- match.arg(family)
  penalty <- match.arg(penalty)

  # error checking: later
  y <- as.matrix(y)
  X_2 <- as.matrix(X_2)
  if(!is.null(X_1)){
    X_1 <- as.matrix(X_1)
  }

  # remove zero columns
  rmid <- which(apply(X_2, 2, sum) == 0)
  if(length(rmid) > 0){
    X_2 <- X_2[, -rmid]
    treemat <- treemat[-rmid, ]
  }

  # tree_guided reparameterization: get all the intermediate quantities
  tree.object <- get_tree_object(X_2, treemat, penalty,
                                 group.weight = group.weight,
                                 feature.weight = feature.weight)
  C_1 <- tree.object$C_1  # for generalized lasso
  CNorm_1 <- tree.object$CNorm_1
  C_2 <- as.matrix(tree.object$C_2)  # for group lasso
  CNorm_2 <- tree.object$CNorm_2
  A <- tree.object$A
  g_idx <- tree.object$g_idx

  # get alpha sequence
  alpha <- modstr$alpha
  nalpha <- length(alpha)
  # get the lambda sequence
  lambda <- modstr$lambda
  nlambda <- modstr$nlambda
  lambda.min.ratio <- modstr$lambda.min.ratio
  if(is.null(lambda)){
    if(is.null(group.weight)){
      pg.max <- max(g_idx[, 3])
      wg <- sqrt(g_idx[, 3]/pg.max)
    }else if(group.weight == 'equal'){
      wg <- rep(1, nrow(g_idx))
    }else{
      wg <- group.weight
    }

    if(is.null(feature.weight)){
      wj <- apply(X_2, 2, FUN = function(x){sqrt(sum(x^2))})
    }else if(feature.weight == 'equal'){
      wj <- rep(1, ncol(X_2))
    }else{
      wj <- feature.weight
    }



    lambda.max <- estlambda(y, X_2, A, g_idx, wg, wj)
    lambda <- exp(seq(log(lambda.max), log(lambda.min.ratio * lambda.max), length.out = nlambda))
  }else{
    nlambda <- length(lambda)
    lambda <- sort(lambda, decreasing = TRUE)
  }

  # get default for initial estimation
  if(is.null(X_1)){
    pall <- ncol(A) + 1
    p1 <- 1
    cov.name <- NULL
  }else{
    pall <- ncol(X_1) + ncol(A) + 1
    p1 <- ncol(X_1) + 1
    cov.name <- colnames(X_1)
  }
  if(is.null(gamma.init)){
    gamma.init <- rep(0, pall)
  }
  # get default for binary weight
  if(family == 'logit'){
    if(is.null(weight))
    weight <- c(1-sum(y)/length(y), sum(y)/length(y))
  }
  # pre-define a matrix to save estimation results
  gammacoef <- array(0, dim = c(pall, nlambda, nalpha))
  groupnorm <- array(0, dim = c(nrow(g_idx), nlambda, nalpha))
  # get estimations at each alpha, lambda combination
  for(i in 1:nalpha){
    alphai <- alpha[i]
    #for(j in 1:nlambda){
    #  lambdaj <- lambda[j]
      gammaij <- TSLA(y, X_1, X_2, C_1, CNorm_1, C_2, CNorm_2, g_idx, A, family,
                      gamma.init, weight, lambda, alphai, control)
      # estimation
      gammacoef[, , i] <- gammaij$gamma.est
      # weighted norm
      gnorms <- apply(gammaij$gamma.est, 2, FUN = function(x){
        cal2norm(C_2 %*% x[(p1 + 1):pall], g_idx, 2)
      })
      groupnorm[, , i] <- gnorms
      #groupnorm[, , i] <- cal2norm(C_2 %*% gammaij$gamma.est[(p1 + 1):pall], g_idx, 2)
    #}

  }
  # output as a list
  bin.name <- colnames(X_2)
  return(list(gammacoef = gammacoef,
              groupnorm = groupnorm,
              lambda.seq = lambda,
              alpha.seq = alpha,
              rmid = rmid, family = family,
              cov.name = cov.name, bin.name = bin.name,
              tree.object = tree.object
              ))

   # the A matrix will not involve zero columns in tree.object
}
