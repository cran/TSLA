# function for tree-guided expansion

### internal function
###===========================================================================
# build extend taxonomy and design matrix for every internal node
getexmat <- function(dorigin){
  # dorgin: design matrix for a tree with height of 2
  p <- ncol(dorigin)
  regrelist <- paste0(colnames(dorigin), collapse='*')
  f <- as.formula(paste('~',regrelist))
  # get design matrix with interaction terms
  dmatfull <- model.matrix(f, data = dorigin)
  # get aggregated variable for this node
  regrelist <- paste0(colnames(dorigin), collapse='|')
  f <- as.formula(paste('~',regrelist))
  dmataggr <- model.matrix(f, data = dorigin)
  # remove intercept and original columns
  namelist <- colnames(dmatfull)[(p+2):ncol(dmatfull)]
  dmatfull <- data.frame(dmatfull[, (p+2):ncol(dmatfull)])
  colnames(dmatfull) <- namelist
  dmataggr <- dmataggr[, 2]
  # calculate order for all generated terms
  order <- 1 + base::lengths(regmatches(colnames(dmatfull), gregexpr(":", colnames(dmatfull))))
  # remove columns which consists of zeros only
  sumzeros <- apply(dmatfull, 2, FUN = sum)
  #n <- nrow(dmatfull)
  nonzeros <- which(sumzeros > 0)
  dmatfull <- dmatfull[, nonzeros]
  order <- order[nonzeros]
  namelist <- namelist[nonzeros]
  # combine everything
  namelist <- c(namelist, colnames(dmataggr))
  order <- c(order, 1)
  exdmat <- as.data.frame(cbind(dmatfull, dmataggr))
  colnames(exdmat) <- namelist
  # the last column is always the aggregate feature
  return(list(exdmat, order))
}



# construct extended taxonomy matrix with required interaction terms
# input: original taxonomy matrix without interactions and original design matrix
#' Tree-guided expansion
#'
#' Give the expanded design matrix and the expanded tree structure
#' by adding interactions in conformity to the structure.
#'
#' This function is used by the TSLA method only when the \code{penalty} is
#' selected as "CL2". The all zero columns produced by the
#' interactions are excluded in the output.
#'
#' For the TSLA method, the signs of the coefficients in the linear
#' constraints depend on the order of the term.
#' To better extend the method in implementation, we apply the signs on
#' the feature vectors instead of the regression coefficients.
#' For example, we use feature vector -\eqn{x_{12}} instead of \eqn{x_{12}}.
#' The expanded design matrix \code{x.expand} from this
#' function is adjusted by the signs. The \code{A} matrix and
#' all the coefficients estimated from the package
#' can be explained correspondingly. We also provide \code{x.expand.adj},
#' \code{A.adj}, and \code{beta.coef.adj}
#' as the quantities with the effects of the signs removed.
#'
#'
#' The input tree structure of the original features needs to be
#' constructed as the following:
#' each row corresponds to a variable at the finest level;
#' each column corresponds to an ordered classification level with the
#' leaf level at the left-most and the root level at the right-most;
#' the entry values in each column are the index of the ancestor node of
#' the variable at that level.
#' As we move from left to right, the number of unique values
#' in the column becomes fewer.
#'
#' @param tmatrix Tree structure of the original features in matrix form.
#' @param dmatrix Original design matrix in matrix form.
#'
#' @return A list.
#' \item{x.expand}{The design matrix after expansion.
#' Each column is multiplied by \eqn{(-1)^{r-1}},
#' where r is the order of the corresponding interaction term.}
#' \item{tree.expand}{The tree structure after expansion.}
#' \item{x.expand.adj}{The design matrix after expansion with the
#' effects of signs removed.}
#' @export
#'
#'
getetmat <- function(tmatrix, dmatrix){
  # read original taxonomy table
  #tmatrix <- read.csv(treefile, header = TRUE)
  # read original design matrix
  #dmatrix <- read.csv(dmatin, header = TRUE)

  # store original predictor names in column order
  originalnames <- colnames(dmatrix)
  colnames(dmatrix) <- paste0(rep('X', ncol(dmatrix)),
                              1, '_', 1:ncol(dmatrix))
  p <- nrow(tmatrix) # number of native leaf nodes
  l <- ncol(tmatrix) # number of tree levels

  tmatrix$type <- rep(1, p) # node type: 0 for internal nodes, 1 for leaf nodes
  tmatrix$id <- colnames(dmatrix) # info about the construction of the variable
  tmatrix$order <- rep(1, p) # used to calculate sign for gamma paramization: 1 or -1
  tmatrix$child <- rep(1, p) # direct child of next level
  #rowid <- p
  #tdimc <- ncol(tmatrix)

  for(i in 2:l){
    t <- length(unique(tmatrix[, i])) # number of internal nodes in level l
    # loop over all internal nodes at the ith level
    for(j in 1:t){
      if(i == 2){
        tmat <- tmatrix[which(tmatrix[, i] == j), ]
        did <- which(tmatrix[, i] == j)
      }else{
        tmat <- tmatrix[which(tmatrix[, i] == j & tmatrix$child == 1), ]
        did <- which(tmatrix[, i] == j & tmatrix$child == 1)
        # only consider direct child
      }
      # singleton at the ith level
      if(nrow(tmat) == 1){
        tmatrix$child[did] <- 1

      }else{
        tmatrix$child[did] <- 0
        dmat <- data.frame(dmatrix[, did])
        colnames(dmat) <- colnames(dmatrix)[did]
        addinter <- getexmat(dorigin = dmat)
        # extend dmatrix
        exdmat <- addinter[[1]]
        nadd <- ncol(exdmat)
        colnames(exdmat)[nadd] <- paste0('X', i, '_', j)
        dmatrix <- cbind(dmatrix, exdmat)
        # extend tmatrix
        levelmat <- data.frame(matrix(rep(as.numeric(tmat[1, 1:l]), nadd),
                                      nrow = nadd, ncol=l, byrow = TRUE))

        for(k in 1:(i-1)){
          levelmat[, k] <- (max(tmatrix[, k])+1):(max(tmatrix[, k])+nadd)
        }
        levelmat[nadd, 1:(i-1)] <- 0
        #levelmat[, 2:i-1] <- 0
        #levelmat[, i-1] <- c((max(tmatrix[, i-1])+1):(max(tmatrix[, i-1])+nadd-1), 0)
        #levelmat[, 1] <- (max(tmatrix[, 1])+1):(max(tmatrix[, 1])+nadd)
        levelmat$type <- c(rep(1, nadd-1), 0)
        levelmat$id <- colnames(exdmat)
        levelmat$order <- addinter[[2]]
        levelmat$child <- c(rep(0, nadd-1), 1)
        colnames(levelmat) <- colnames(tmatrix)
        tmatrix <- rbind(tmatrix, levelmat)


      }
    }
  }

  did <- which(tmatrix[, 1] != 0)
  finaltmat <- tmatrix[tmatrix$type == 1, ]
  finaldmat <- dmatrix[, did]
  # add sign to the design matrix
  for(i in 1:ncol(finaldmat)){
    finaldmat[, i] <- finaldmat[, i]*(-1)^(finaltmat$order[i]-1)
  }
  # add label to nonzero predictors(need to adjust based on the tree mannually)
  #finaltmat$nonzero <- 0
  #finaltmat$
  #nonzero[which(finaltmat$Order==1|finaltmat$Order==3|finaltmat$Class==2)] <- 1


  #write.table(finaltmat, tmatout, col.names = TRUE,
  #            row.names = FALSE, sep = ',')
  #write.table(finaldmat, dmatout, col.names = TRUE,
  #            row.names = FALSE, sep = ',')
  # remove intermediate columns
  finaltmat <- finaltmat[, 1:l]
  return(list(x.expand = as.matrix(finaldmat),
              tree.expand = as.matrix(finaltmat),
              x.expand.adj = as.matrix(abs(finaldmat))))
}




