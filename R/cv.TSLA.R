#source("./TSLA.fit.R")
#source("./tools.R")
# test
#' Cross validation for TSLA
#'
#' Conduct cross validation to select the optimal tuning parameters in TSLA.
#'
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
#' @param pred.loss Model performance metrics. If \code{family="ls"}, default
#' is "MSE" (mean squared error). If \code{family="logit"}, default is "AUC". For logistic
#' model, another option is "deviance".
#' @param gamma.init Initial value for the optimization. Default is a zero vector.
#' The length should equal to 1+\code{ncol(X_1)}+\code{ncol(A)}.
#' See details of A in \code{get_tree_obj()}.
#' @param weight A vector of length two and it is used for logistic regression
#' only. The first element corresponds to weight of y=1 and the
#' second element corresponds to weight of y=0.
#' @param nfolds Number of cross validation folds. Default is 5.
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
#' @return A list of cross validation results.
#' \item{lambda.min}{\eqn{\lambda} value with best prediction performance.}
#' \item{alpha.min}{\eqn{\alpha} value with best prediction performance.}
#' \item{cvm}{A (number-of-lambda * number-of-alpha) matrix saving the means of cross validation loss across folds.}
#' \item{cvsd}{A (number-of-lambda * number-of-alpha) matrix saving standard deviations of
#'  cross validation loss across folds.}
#'  \item{TSLA.fit}{Outputs from \code{TSLA.fit()}.}
#'  \item{Intercept.min}{Intercept corresponding to \code{(lambda.min,alpha.min)}.}
#'  \item{cov.min}{Coefficients of unpenalized features
#'  corresponding to \code{(lambda.min,alpha.min)}.}
#'  \item{beta.min}{Coefficients of binary features corresponding
#'  to \code{(lambda.min,alpha.min)}.}
#'  \item{gamma.min}{Node coefficients corresponding to \code{(lambda.min,alpha.min)}.}
#'  \item{groupnorm.min}{Group norms of node coefficients corresponding to \code{(lambda.min,alpha.min)}.}
#'  \item{lambda.min.index}{Index of the best \eqn{\lambda} in the sequence.}
#'  \item{alpha.min.index}{Index of the best \eqn{\alpha} in the sequence.}
#'
#' @export
#'
#'
#' @examples
#' # Load the synthetic data
#' data(ClassificationExample)
#'
#' tree.org <- ClassificationExample$tree.org   # original tree structure
#' x2.org <- ClassificationExample$x.org      # original design matrix
#' x1 <- ClassificationExample$x1
#' y <- ClassificationExample$y            # response
#'
#' # Do the tree-guided expansion
#' expand.data <- getetmat(tree.org, x2.org)
#' x2 <- expand.data$x.expand              # expanded design matrix
#' tree.expand <- expand.data$tree.expand  # expanded tree structure
#'
#' # Do train-test split
#' idtrain <- 1:200
#' x1.train <- as.matrix(x1[idtrain, ])
#' x2.train <- x2[idtrain, ]
#' y.train <- y[idtrain, ]
#' x1.test <- as.matrix(x1[-idtrain, ])
#' x2.test <- x2[-idtrain, ]
#' y.test <- y[-idtrain, ]
#'
#' # specify some model parameters
#' set.seed(100)
#' control <- list(maxit = 100, mu = 1e-3, tol = 1e-5, verbose = FALSE)
#' modstr <- list(nlambda = 5,  alpha = seq(0, 1, length.out = 5))
#' simu.cv <- cv.TSLA(y = y.train, as.matrix(x1[idtrain, ]),
#'                    X_2 = x2.train,
#'                    treemat = tree.expand, family = 'logit',
#'                    penalty = 'CL2', pred.loss = 'AUC',
#'                    gamma.init = NULL, weight = c(1, 1), nfolds = 5,
#'                    group.weight = NULL, feature.weight = NULL,
#'                    control = control, modstr =  modstr)
#' # Do prediction with the selected tuning parameters on the test set. Report AUC on the test set.
#' rmid <- simu.cv$TSLA.fit$rmid  # remove all zero columns
#' if(length(rmid) > 0){
#'   x2.test <- x2.test[, -rmid]}
#'   y.new <- predict_cvTSLA(simu.cv, as.matrix(x1[-idtrain, ]), x2.test)
#'   library(pROC)
#'   auc(as.vector(y.test), as.vector(y.new))
#'
#'
cv.TSLA <- function(y, X_1 = NULL, X_2, treemat, family = c('ls', 'logit'),
                    penalty = c('CL2', 'RFS-Sum'),
                    pred.loss = c('MSE', 'AUC', 'deviance'),
                    gamma.init = NULL, weight = NULL, nfolds = 5,
                    group.weight = NULL, feature.weight = NULL,
                    control = list(), modstr = list()){
  ## X_1: covariate matrix
  ## X_2: expanded design matrix
  Call <- match.call()
  ## set control values for model fitting
  control <- do.call("TSLA.control", control)
  ## set mod strings
  modstr <- do.call("TSLA.modstr", modstr)
  ## get choice of penalty form
  penalty <- match.arg(penalty)
  ## make sure default for logistic regression is AUC
  if(family == 'logit' & pred.loss == 'MSE'){
    pred.loss <- 'AUC'
  }
  ## get choice of loss function
  pred.loss <- match.arg(pred.loss)

  y <- as.matrix(y)
  X_2 <- as.matrix(X_2)
  if(!is.null(X_1)){
    X_1 <- as.matrix(X_1)
  }

  # run TSLA.fit for the first time to get the whole solution path
  TSLA.object <- TSLA.fit(y, X_1 = X_1, X_2 = X_2, treemat = treemat,
                          family = family, penalty = penalty,
                          gamma.init = gamma.init, weight = weight,
                          group.weight = group.weight,
                          feature.weight = feature.weight,
                          control = control, modstr = modstr)
  ## get lambda and alpha range
  # get alpha sequence
  alpha <- TSLA.object$alpha.seq
  nalpha <- length(alpha)
  # get the lambda sequence
  lambda <- TSLA.object$lambda.seq
  nlambda <- modstr$nlambda
  modstr$lambda <- lambda
  modstr$alpha <- alpha
  n <- nrow(X_2)

  # do cross-validation
  if(family == 'logit'){
    foldid0 <- sample(rep(seq(nfolds), length = n - sum(y)))
    foldid1 <- sample(rep(seq(nfolds), length = sum(y)))
    foldid <- rep(0, n)
    foldid[which(y == 0)] <- foldid0
    foldid[which(y == 1)] <- foldid1
  }else{
    foldid <- sample(rep(seq(nfolds), length = n))
  }

  outlist <- as.list(seq(nfolds))
  for(i in seq(nfolds)){
    # train/test splitting
    #if(family == 'logit'){
    #  trainid <- c(which(foldid0 != i), which(foldid1 != i))
    #}else{
      trainid <- which(foldid != i)
    #}

    if(!is.null(X_1)){
      X_1.train <- X_1[trainid, ]
    }else{
      X_1.train <- NULL
    }

    y.train <- as.matrix(y)[trainid, ]
    X_2.train <- X_2[trainid, ]
    outlist[[i]] <- TSLA.fit(y.train, X_1.train, X_2.train, treemat,
                            family = family, penalty = penalty,
                         gamma.init = gamma.init, weight = weight,
                         group.weight = group.weight,
                         feature.weight = feature.weight,
                         control = control, modstr = modstr)
  }
  ###What to do depends on the pred.loss and family
  fun <- paste("cv", family, sep = ".")
  cvstuff <- do.call(fun, list(outlist, y, X_1, X_2, foldid, pred.loss))
  cvm <- cvstuff$cvm
  cvsd <- cvstuff$cvsd
  paramin <- cvstuff$paramin
  if(nrow(paramin) == 1){
    lambda.min.index <- paramin[1, 1]
    alpha.min.index <- paramin[1, 2]
    lambda.min <- lambda[lambda.min.index]
    alpha.min <- alpha[alpha.min.index]
  }else{
    lambda.min.index <- min(paramin[, 1])
    col.index <- which(paramin[, 1] == min(paramin[, 1]))
    lambda.min <- lambda[min(paramin[, 1])]
    # max alpha to get more features aggregated
    alpha.min.index <- paramin[, 2][col.index ]
    alpha.min.index  <- max(alpha.min.index)
    alpha.min <- alpha[alpha.min.index ]
  }
  ### get corresponding coefficient from TSLA.fit
  coefs <- coef_TSLA(TSLA.object)
  Intercept.min <- coefs$Intercept[, lambda.min.index, alpha.min.index]
  if(is.null(X_1)){
    cov.min <- NULL
  }else if(ncol(X_1) == 1){
    cov.min <- coefs$cov.coef[, lambda.min.index, alpha.min.index]
  }else{
    cov.min <- coefs$cov.coef[, lambda.min.index, alpha.min.index]
  }
  beta.min <- coefs$beta.coef[, lambda.min.index, alpha.min.index]
  gamma.min <- coefs$gamma.coef[, lambda.min.index, alpha.min.index]
  groupnorm.min <- TSLA.object$groupnorm[, lambda.min.index, alpha.min.index]
  out <- list(lambda.min = lambda.min, alpha.min = alpha.min,
              cvm = cvm, cvsd = cvsd,
              TSLA.fit = TSLA.object,
              Intercept.min = Intercept.min,
              cov.min = cov.min, beta.min = beta.min,
              gamma.min = gamma.min,
              groupnorm.min = groupnorm.min,
              lambda.min.index = lambda.min.index,
              alpha.min.index = alpha.min.index)
  return(out)
}

### internal functions
# sub-functions to calculate loss for each fold
# default for ls: MSE
cv.ls <- function(outlist, y, X_1 = NULL, X_2, foldid, pred.loss = 'MSE'){
  pred.loss <- match.arg(pred.loss)
  nfolds <- max(foldid)
  fitdim <- dim(outlist[[1]]$gammacoef)
  out <- array(0, dim = c(nfolds, fitdim[2], fitdim[3]))
  for(i in seq(nfolds)){
    fitobj <- outlist[[i]]
    testid <- which(foldid == i)
    ni <- length(testid)
    y.test <- as.matrix(y)[testid, ]
    if(!is.null(X_1)){
      X_1.test <- as.matrix(X_1[testid, ])
    }else{
      X_1.test <- NULL
    }
    # remove zero columns identified in the training set
    X_2.test <- X_2[testid, ]
    if(length(fitobj$rmid) > 0){
      X_2.test <- X_2.test[, -fitobj$rmid]
    }
    # get prediction at each lambda, alpha combination in the ith fold
    ypre <- predict_TSLA(fitobj, X_1.test, X_2.test, type = "response")
    out[i, , ] <- apply(ypre, c(2, 3), FUN = function(x) sum((x-y.test)^2))/ni
  }
  cvm <- apply(out, c(2, 3), FUN = mean)
  cvsd <- apply(out, c(2, 3), FUN = sd)
  paramin <- which(cvm == min(cvm), arr.ind = TRUE)
  return(list(cvraw = out, cvm = cvm, cvsd = cvsd, paramin = paramin))
}

# sub-functions to calculate loss for each fold
# default for logistic: AUC
cv.logit <- function(outlist, y, X_1 = NULL, X_2, foldid,
                     pred.loss = c('AUC', 'deviance')){
  pred.loss <- match.arg(pred.loss)
  nfolds <- max(foldid)
  fitdim <- dim(outlist[[1]]$gammacoef)
  out <- array(0, dim = c(nfolds, fitdim[2], fitdim[3]))
  for(i in seq(nfolds)){
    fitobj <- outlist[[i]]
    testid <- which(foldid == i)
    ni <- length(testid)
    y.test <- as.matrix(y)[testid, ]
    if(!is.null(X_1)){
      X_1.test <- as.matrix(X_1[testid, ])
    }else{
      X_1.test <- NULL
    }
    # remove zero columns identified in the training set
    X_2.test <- X_2[testid, ]
    if(length(fitobj$rmid) > 0){
      X_2.test <- X_2.test[, -fitobj$rmid]
    }
    # get prediction at each lambda, alpha combination in the ith fold
    ypre <- predict_TSLA(fitobj, X_1.test, X_2.test, type = "response")
    outi <- switch(pred.loss,
                   "AUC" = apply(ypre, c(2, 3), FUN = function(x) as.numeric(auc(y.test, x))),
                   "deviance" = apply(ypre, c(2, 3), FUN = function(x){
                     sum(-log(x[y.test == 1])) + sum(-log(1 - x[y.test == 0]))
                   } ))
    out[i, , ] <- outi
  }
  cvm <- apply(out, c(2, 3), FUN = mean)
  cvsd <- apply(out, c(2, 3), FUN = sd)
  paramin <- switch(pred.loss,
                    "AUC" = which(cvm == max(cvm), arr.ind = TRUE),
                    "deviance" = which(cvm == min(cvm), arr.ind = TRUE))

  return(list(cvraw = out, cvm = cvm, cvsd = cvsd, paramin = paramin))
}
