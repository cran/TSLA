##library(pROC)
#library(PRROC)
## function that generate coefficients
## modified from gglasso package
## get coefficient from a TSLA.fit object
#' Get coefficients from a fitted TSLA model
#'
#' Get coefficients from a TSLA.fit object.
#' @param object A fit output from \code{TSLA.fit()}.
#'
#' @param ... Other parameters.
#' @return A list of coefficients for each combination of \eqn{\lambda} and
#' \eqn{\alpha}.
#' The first dimension is indexed by the coefficient,
#' the second dimension is indexed by \eqn{\lambda},
#' and the third dimension is
#' indexed by \eqn{\alpha}.
#' \item{Intercept}{Intercept.}
#' \item{cov.coef}{Coefficients for unpenalized features.}
#' \item{gamma.coef}{Node coefficients. Also see details in \code{getetmat()}.}
#' \item{beta.coef}{Regression coefficients. Each coefficient is multiplied by
#' \eqn{(-1)^{r-1}}, where r is the order of the corresponding interaction term.
#' Also see details in \code{getetmat()}.}
#' \item{beta.coef.adj}{Regression coefficients \eqn{\beta}.
#' This is the common regression coefficient vector.
#' It corresponds to \code{x.expand.adj} and
#' \code{A.adj}.}

#' @export
#'
coef_TSLA <- function(object, ...){
  TSLA.object <- object
  A <- TSLA.object$tree.object$A
  n.beta <- nrow(A)
  n.gamma <- ncol(A)
  p <- dim(TSLA.object$gammacoef)[1]
  n.cov <-  p - n.gamma - 1
  n.lambda <- length(TSLA.object$lambda.seq)
  n.alpha <- length(TSLA.object$alpha.seq)
  # intercept
  intercept <- array(TSLA.object$gammacoef[1, , ], c(1, n.lambda, n.alpha))
  # unpenalized covariates
  if(n.cov == 0){
    cov.coef <- NULL
  }else if(n.cov == 1){
    cov.coef <- array(TSLA.object$gammacoef[2, , ], c(1, n.lambda, n.alpha))
  }else{
    if((n.cov > 1) & is.null(TSLA.object$cov.name)){
      cov.name <- paste0('cov.', 1:n.cov)
      cov.coef <- array(TSLA.object$gammacoef[2:(n.cov + 1), , ], c(n.cov, n.lambda, n.alpha))
      dimnames(cov.coef)[[1]] <- cov.name
    }else{
      cov.name <- TSLA.object$cov.name
      cov.coef <- array(TSLA.object$gammacoef[2:(n.cov + 1), , ], c(n.cov, n.lambda, n.alpha))
      dimnames(cov.coef)[[1]] <- cov.name
    }
  }

  # gamma coefficient
  # for gamma coefficient, use node id as name, at least 2 nodes
  node.name <- paste0('node.', 1:n.gamma)
  gamma.coef <- array(TSLA.object$gammacoef[(n.cov+2):p, , ], c(n.gamma, n.lambda, n.alpha))
  dimnames(gamma.coef)[[1]] <- node.name

  # for binary features: at least 2
  if(is.null(TSLA.object$bin.name)){
    bin.name <- paste0('binary.', 1:n.beta)
  }else{
    bin.name <- TSLA.object$bin.name
  }
  beta.coef <- apply(gamma.coef, c(2, 3), FUN = function(x) A %*% x)
  dimnames(beta.coef)[[1]] <- bin.name

  # for binary features: at least 2
  beta.coef.adj <- apply(gamma.coef, c(2, 3), FUN = function(x) TSLA.object$tree.object$A.adj %*% x)
  dimnames(beta.coef.adj)[[1]] <- bin.name

  ## when s is a given combination of lambda and alpha?


  return(list(Intercept = intercept,
              cov.coef = cov.coef,
              gamma.coef = gamma.coef,
              beta.coef = beta.coef,
              beta.coef.adj = beta.coef.adj))
}


## function to do prediction with new data points, also for crossvalidation evaluation
#' Prediction from TSLA with new data
#'
#' Generate prediction for the response.
#' @param object A fit output from \code{TSLA.fit()}.
#' @param X_1_new New unpenalized features in matrix form.
#' @param X_2_new New binary features in matrix form.
#' @param type Two options: "response" or "link". The two options only
#' differ for \code{family="logit"}.
#' @param ... Other parameters.
#'
#' @return Predictions. The first dimension is indexed by the observation,
#' the second dimension is indexed by \eqn{\lambda},
#' and the third dimension is indexed by \eqn{\alpha}.

#' @export
predict_TSLA <- function(object, X_1_new = NULL, X_2_new,
                         type = c('response', 'link'), ...){
  TSLA.object <- object
  type <- match.arg(type)
  coefs <- coef_TSLA(TSLA.object)
  # add intercept
  allone <- matrix(rep(1, nrow(X_2_new)), ncol = 1)
  if(is.null(X_1_new)){
    nfit <- apply(coefs$Intercept, c(2, 3), FUN = function(x) allone %*% x) +
            apply(coefs$beta.coef, c(2, 3), FUN = function(x) X_2_new %*% x)
  }else if(ncol(X_1_new) == 1){
    nfit <- apply(coefs$Intercept, c(2, 3), FUN = function(x) allone %*% x) +
      apply(coefs$cov.coef, c(2, 3), FUN = function(x) X_1_new %*% x) +
      apply(coefs$beta.coef, c(2, 3), FUN = function(x) X_2_new %*% x)
  }else{
    nfit <- apply(coefs$Intercept, c(2, 3), FUN = function(x) allone %*% x) +
            apply(coefs$cov.coef, c(2, 3), FUN = function(x) X_1_new %*% x) +
            apply(coefs$beta.coef, c(2, 3), FUN = function(x) X_2_new %*% x)
  }
  family <- TSLA.object$family
  if(family == 'ls'){
    return(nfit)
  }else{
    nfit <- switch(type,
           'link' = nfit,
           'response' = apply(nfit, c(2, 3), FUN = function(x) exp(x)/(1+exp(x))))
    return(nfit)
  }
}


## function to do prediction from cross-validation results
#' Prediction from cross validation
#'
#' A convenient function to get prediction from the
#' selected tuning parameters by cross validation.
#'
#' @param object A fit output from \code{cv.TSLA().}
#'
#' @param X_1_new New unpenalized features in matrix form.
#' @param X_2_new New binary features in matrix form.
#' @param type Two options: "response" or "link". The two options only
#' differ for \code{family="logit"}.
#' @param ... Other parameters.
#'
#' @return Predictions.
#' @export
predict_cvTSLA <- function(object, X_1_new = NULL, X_2_new,
                         type = c('response', 'link'), ...){
  cv.TSLA.object <- object
  type <- match.arg(type)
  Intercept <- cv.TSLA.object$Intercept.min
  beta <- cv.TSLA.object$beta.min
  cov <- cv.TSLA.object$cov.min
  # add intercept
  allone <- matrix(rep(1, nrow(X_2_new)), ncol = 1)
  if(is.null(X_1_new)){
    nfit <- allone * Intercept + X_2_new %*% beta
  }else if(ncol(X_1_new) == 1){
    nfit <- allone * Intercept + X_1_new * cov + X_2_new %*% beta
  }else{
    nfit <- allone * Intercept + X_1_new %*% cov + X_2_new %*% beta
  }
  family <- cv.TSLA.object$TSLA.fit$family
  if(family == 'ls'){
    return(nfit)
  }else{
    nfit <- switch(type,
                   'link' = nfit,
                   'response' = exp(nfit)/(1+exp(nfit)))
    return(nfit)
  }
}

## function to evaluate model performance

#' Get performance metrics for classification
#'
#' Evaluate the prediction performance under the classification settings.
#'
#' The function supports three methods to select the threshold of the
#' predicted probability.
#'
#' \code{threshold.method = "youden"}: The optimal threshold corresponds to
#' the point that maximizes the distance to the identity (diagonal) line on
#' the ROC curve.
#'
#' \code{threshold.method = "specificity.control"}: The optimal threshold
#' corresponds to the smallest value that ensures the required specificity
#' value.
#'
#' \code{threshold.method = "quantile"}: The optimal threshold corresponds to
#' the required quantile of the predicted probability.
#'
#'
#' @param ytest Response vector for test data.
#' @param ypretest Predicted probability for test data.
#' @param family "ls" or "logic". Return MSE when "ls" is used.
#' @param threshold.method Method to get the threshold.
#' @param specificity User-defined specificity or quantile.
#'
#' @return List of measures.
#' \item{AUC}{Area under the ROC curve.}
#' \item{AUPRC}{Area under the precision-recall curve.}
#' \item{threshold}{Selected threshold of the probability.}
#' \item{sensitivity}{Sensitivity with the selected threshold.}
#' \item{ppv}{Positive predictive value with the selected threshold.}
#' \item{specificity}{Specificity with the selected threshold.}
#' \item{true.positive}{Number of true positive with the selected threshold.}
#' \item{false.positive}{Number of false positive with the selected threshold.}
#' @export
getperform <- function(ytest, ypretest, family,
                       threshold.method =
                         c('youden', 'specificity.control', 'quantile'),
                       specificity = NULL){
  if(family == 'ls'){
    model.loss <- sum((ytest - ypretest)^2)/nrow(as.matrix(ytest))
    return(list(MSE = model.loss))
  }
  if(family == 'logit'){
    # get AUC and AUPRC
    AUC <- auc(ytest, ypretest)
    pr <- pr.curve(scores.class0 = ypretest[ytest == 1],
                   scores.class1 = ypretest[ytest == 0],
                   curve=T)
    AUPRC <- pr$auc.integral
    # get sensitivity, PPV, specificity
    threshold.method <- match.arg(threshold.method)
    if(threshold.method == 'youden'){
      roccurve <- roc(ytest, ypretest)
      threshold <- as.numeric(min(coords(roccurve, "best", ret = "threshold")))
      prey <- rep(0, length(ypretest))
      prey[which(ypretest > threshold)] <- 1
      sensitivity <- sum(ytest == 1 & prey == 1)/sum(ytest)
      ppv <- sum(ytest == 1 & prey == 1)/sum(prey)
      specificity <- sum(ytest == 0 & prey ==0)/sum(ytest == 0)
      tp <- sum(ytest == 1 & prey == 1)
      fp <- sum(ytest == 0 & prey == 1)
    }else if(threshold.method == 'specificity.control'){
      roccurve <- roc(ytest, ypretest)
      rocresult <- data.frame(sens = roccurve$sensitivities,
                              spec = roccurve$specificities,
                              thre = roccurve$thresholds)
      # find threshold most close to specified specificity value
      idt <- which(rocresult$spec > specificity)
      rocresult <- rocresult[idt, ]
      idt <- which(abs(rocresult$spec-specificity) == min(abs(rocresult$spec-specificity)))
      rocresult <- rocresult[idt, ]
      # if multiple point, select the one closest to upper left
      idt <- which(rocresult$spec+rocresult$sens == max(rocresult$spec+rocresult$sens))
      rocresult <- rocresult[idt, ]
      # if still multiple, select smallest threshold
      threshold <- min(rocresult$thre)
      prey <- rep(0, length(ypretest))
      prey[which(ypretest > threshold)] <- 1
      sensitivity <- sum(ytest == 1 & prey == 1)/sum(ytest)
      ppv <- sum(ytest == 1 & prey == 1)/sum(prey)
      specificity <- sum(ytest == 0 & prey ==0)/sum(ytest == 0)
      tp <- sum(ytest == 1 & prey == 1)
      fp <- sum(ytest == 0 & prey == 1)
    } else {
      threshold <- quantile(ypretest, specificity)
      prey <- rep(0, length(ypretest))
      prey[which(ypretest > threshold)] <- 1
      sensitivity <- sum(ytest == 1 & prey == 1)/sum(ytest)
      ppv <- sum(ytest == 1 & prey == 1)/sum(prey)
      specificity <- sum(ytest == 0 & prey ==0)/sum(ytest == 0)
      tp <- sum(ytest == 1 & prey == 1)
      fp <- sum(ytest == 0 & prey == 1)
    }
    return(list(AUC = AUC, AUPRC = AUPRC, threshold = threshold,
                sensitivity = sensitivity, ppv = ppv, specificity = specificity,
                true.positive = tp, false.postive = fp))


  }
}
