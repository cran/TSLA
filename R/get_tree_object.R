## function used to get model values for SPG modeling
#source('./reparameterization.R')
#library(Matrix)
#' Tree-guided reparameterization
#'
#' This function generates all the intermediate quantities
#' based on the tree-guided reparameterization.
#'
#' @param X_2 Expanded design matrix for \code{penalty = "CL2"};
#' Original design matrix
#' for \code{penalty = "RFS-Sum"}. Need to be in the matrix form.
#' @param treemat Expanded tree structure for
#' \code{penalty = "CL2"}; Original structure for \code{penalty = "RFS-Sum"}.
#' Need to be in the matrix form.
#' @param penalty Two options for group penalty on \eqn{\gamma},
#' "CL2" or "RFS-Sum".
#' @param group.weight User-defined weights for group penalty.
#' Need to be a vector and the length equals to the number of groups.
#' @param feature.weight User-defined weights for each predictor
#' after expansion.
#'
#' @return A list consists of quantities needed for SPG optimization.
#' \item{C_1}{C_1 matrix for generalized lasso penalty.}
#' \item{CNorm_1}{Nuclear norm of matrix \code{C_1}.}
#' \item{C_2}{C_2 matrix for group lasso penalty.}
#' \item{CNorm_2}{Nuclear norm of matrix \code{C_2}.}
#' \item{A}{A (number-of-leaf * number-of-node) binary matrix containing linear constraints.
#' Recall that \eqn{\beta=A\gamma}.
#' It is used with \code{beta.coef} and \code{x.expand}.}
#' \item{g_idx}{A (number-of-group * 3) matrix.
#' Each column stands for starting row in \code{C_2} of a group,
#' end row in \code{C_2} of a group, and the group size.}
#' \item{M2}{A (number-of-leaf * number-of-level) node index matrix,
#' with index going from 1 to the number of nodes.
#' Root node has index equal to the number of nodes.
#' Each row corresponds to a variable at the finest level,
#' each column corresponds to an ordered classification level;
#' the entry values in each column are the unique indices of the variables
#' at that level.
#' As we move to the right, the number of unique values becomes fewer.}
#' \item{Tree}{A (number-of-group *  number-of-node) group index matrix.
#' Each row is a group and the column order is the same as the order of node
#' index in M2. If the jth node belongs to the ith group, then the (i, j)
#' element of the matrix is 1; otherwise the element is 0.}
#' \item{A.adj}{A (number-of-leaf * number-of-node) binary matrix
#' containing linear constraints.
#' It is used with \code{beta.coef.adj} and \code{x.expand.adj}.}
#'
#' @export
#'
get_tree_object <- function(X_2, treemat, penalty = c('CL2', 'DL2', 'RFS-Sum'),
                            group.weight = NULL,
                            feature.weight = NULL){
  penalty <- match.arg(penalty)
  # tree-guided reparameterization
  group_object <- mat2SPGtree(treemat, penalty, group.weight)

  # for group lasso
  if(penalty == 'RFS-Sum'){
    A <- group_object$A
    if(is.null(group.weight)){
      X_gamma <- X_2 %*% A
      wl <- apply(X_gamma, 2, FUN=function(x){sqrt(sum(x^2))})
      n <- nrow(X_gamma)
      D <- diag(wl)/sqrt(n)
    }else if(group.weight == 'equal'){
      D <- diag(c(rep(1, ncol(A)-1), 0))
    }else{
      D <- diag(group.weight)
    }
    C_2 <- D * group_object$C
    CNorm_2 <- max(colSums(as.matrix(C_2^2)));
  }else{
    C_2 <-group_object$C
    CNorm_2 <- group_object$CNorm
  }

  # for generalized lasso
  if(is.null(feature.weight)){
    wj <- apply(X_2, 2, FUN=function(x){sqrt(sum(x^2))})
    n <- nrow(X_2)
    D <- diag(wj)/sqrt(n)
  }else if(feature.weight == 'equal'){
    D <- diag(rep(1, ncol(X_2)))
  }else{
    D <- diag(feature.weight)
  }


  C_1 <- D %*% group_object$A
  CNorm_1 <- max(colSums(C_1^2))
  #CNorm_1 <- norm(Matrix(C_1), type = "2")^2

  # signs
  colcounts <- apply(X_2, 2, sum)
  signs <- sign(colcounts)
  A.adj <- diag(signs) %*% group_object$A

  return(list(C_1 = C_1, CNorm_1 = CNorm_1, C_2 = as.matrix(C_2),
              CNorm_2 = CNorm_2,
              A = group_object$A, g_idx = group_object$g_idx, M2  = group_object$M2,
              Tree = group_object$Tree, A.adj = A.adj))
}


