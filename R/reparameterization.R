## internal functions
## This function translate a tree structure (parent node vector without redundant nodes from mat2node function)
## to two sets of parameters: the group membership parameter C, CNorm, g.idx (induced from T and Tw), and tree guided model expansion matrix A.
## input:
##    M: p * L matrix, each row corresponds to a variable at the
##       finest level, each column corresponds to an ordered classification level; the
##       entry values in each column are the unique ID of the variable
##       at that level. As we move to the right, the # of unique values become fewer.
##   penalty: the indicator for the form of penalty function.
##            1 (default): overlapping group penalty, L2 of all descendant nodes for each internal node;
##            2: non-overlapping group penalty, L2 of direct child nodes for each internal node;
##            3: lasso penalty for each node (excluding root node)
## ouput:
##   A: #leaf * (#node) binary expansion matrix (for tree regression reparameterization)
##   C: sum(group size) * (#node), very tall, sparse matrix (for SPG)
##   CNorm: the norm of C, as defined in the SPG paper (for SPG)
##   g.idx: #groups * 3 matrix, starting row in C of a group, end row of a group, group size
##   node: 1 * #node vector, [output of mat2node], a parent node vector with p leaf nodes; NO redundant internal
##         nodes; use treeplot to show the actural tree structure
##   M2: p * L or p * (L+1) node index matrix, [output of mat2node], similar to M, but cleaner,
##       with index going from 1 to #nodes (root node has index equal to #node).
##   Tree: #groups * #node matrix, group index matrix, each row is
##      a group, the column order is the same as in node or M2
##      For instance, if penalty=1 or 2, T has size (#node-#leaf) * #node;
##      if penalty=3, T has size #node-1 * #node
##   Tw: length(nrow(T)) vector,  weight for each group, default is a vector of 1
mat2SPGtree <- function(M, penalty = c("CL2", "DL2", "RFS-Sum"), group.weight = NULL) {
  penalty <- match.arg(penalty)
  ## call mat2node to get tree vector
  tmp.node <- mat2node(M)
  node <- tmp.node$node
  M2 <- tmp.node$M2

  p <- dim(M)[1]  ## number of coefficient
  numall <- length(node)  ## number of nodes in the tree T
  numint <- numall - p   ## number of internal nodes (including root)

  ## expansion matrix A from tree-guided reparameterization
  A <- matrix(0, nrow = p, ncol = numall)
  for (j in 1:p)  { ## leaf node/variable in model
    ## first, get all parents of the leaf node
    parent <- unlist(unique(M2[j, ]))
    A[j, parent] <- 1
  }

  ## get the membership matrix T1
  if(penalty == "DL2") {  ## overlapping group lasso of all descendant nodes
    T1 <- matrix(0, nrow = numint + 1, ncol = numall)
    T2 <- matrix(0, nrow = numint + 1, ncol = numall)
    for(j in 1:numint) {
      nodeindex <- j + p  ## index of internal nodes
      ixy <- which(M2 == nodeindex, arr.ind=TRUE)
      iy1 <- max(ixy[, 2])
      ix1 <- ixy[ixy[, 2] == iy1, 1]
      child <- setdiff(M2[ix1, 1:iy1], nodeindex) ## all descendants of nodeindex
      T1[j, child] <- 1
      iy1 <- min(ixy[, 2])
      ix1 <- ixy[ixy[, 2] == iy1, 1]
      child <- unlist(unique(M2[ix1, iy1-1])) ## direct child of the nodeindex
      T2[j, child] <- 1
    }
    T1[numint + 1, ] <- 1  # add root node as a singleton
    T2[numint + 1, ] <- 1
  } else if(penalty == "CL2") {  ## group lasso of direct child nodes
    T1 <- matrix(0, nrow = numint + 1, ncol = numall)
    for(j in 1:numint) {
      nodeindex <- j + p  ## index of internal nodes
      ixy <- which(M2 == nodeindex, arr.ind=TRUE)
      iy1 <- min(ixy[, 2])
      ix1 <- ixy[ixy[, 2] == iy1, 1]
      child <- unlist(unique(M2[ix1, iy1 - 1])) ## all descendants of nodeindex
      T1[j, child] <- 1
    }
    T1[numint + 1, numall] <- 1
    T2 <- T1
  } else {
    T1 <- diag(numall)
    T2 <- T1
  }

  ## weight for different node groups
  if(is.null(group.weight)){
    k <- max(apply(T2, 1, sum))
    group.weight <- switch(penalty,
                          "CL2" = sqrt(apply(T2, 1, sum)/k),
                          "DL2" = sqrt(apply(T2, 1, sum)/k),
                          "RFS-Sum" = rep(1, nrow(T2)))
  }else if(group.weight == 'equal'){
    group.weight <- rep(1, nrow(T2))
  }else{
    group.weight <- group.weight
  }
  #group.weight <- rep(1, nrow(T2))



  ## SPG required input
  C.res <- pre_group(T1, group.weight)
  return(list(
    A = A, C = C.res$C, CNorm = C.res$CNorm, g_idx = C.res$g_idx,
    node = node, M2 = M2, Tree = T1, group.weight = group.weight))

}


## the pre grouping function from SPG
pre_group <- function(T1, group.weight) {
  V <- dim(T1)[1]
  K <- dim(T1)[2]
  sum_col_T <- rowSums(T1)
  SV <- sum(sum_col_T)
  csum <- cumsum(sum_col_T)
  g_idx <- cbind(c(1, head(csum, -1) + 1), csum, sum_col_T) ## each row is the range of the group
  J <- rep(0, SV)
  W <- rep(0, SV)
  for(v in 1:V) {
    J[g_idx[v, 1]:g_idx[v, 2]] <- which(T1[v, ] == 1)
    W[g_idx[v, 1]:g_idx[v, 2]] <- group.weight[v];
  }
  C <- Matrix::sparseMatrix(i = 1:SV, j = J, x = W, dims = c(SV, K))
  TauNorm <- max(colSums((diag(group.weight) %*% T1)^2))
  return(list(C = C, CNorm = TauNorm, g_idx = g_idx))
}


# this function converts a matrix M representing tree structure to a tree
# parent-node vector and a standarded matrix M2 (with node index from 1 to numnode)
mat2node <- function(M) {
  p <- dim(M)[1] ## number of leaf node
  l <- dim(M)[2] ## number of classification levels
  ## check if the finest level has unique taxa
  if(length(unique(M[, 1])) < p) {
    stop('The first level has overlapping taxa! Terminated.');
  }
  ## check if M is really a tree
  for(j in 1:(l-1)) {
    uniqid <- unique(M[,j])
    for(k in 1:length(uniqid)) {
      if (length(unique(M[M[,j]==uniqid[k],j+1]))>1) { ## check if nested
        stop('The input matrix does not have a tree structure! Terminated.')
      }
    }
  }
  ## if the last level is not all equal, add a root node
  if (length(unique(M[, l]))!=1) {
    # warning('Adding a root node on top of the highest classification level.')
    M <- cbind(M, 1)
    l <- l + 1
  }
  ## Parent node vector with redundant nodes in the tree
  M1 <- M ## node index matrix, index goes from 1 to # of nodes
  M1[, 1] <- 1:p
  for(j in 2:l) {
    M1[, j] <- match(M[, j], unique(M[, j]), nomatch = 0)
    ## now the numbers are unique node indices, the last column values should be the total number of nodes in the tree
    M1[, j] <- M1[, j] + max(M1[, j-1])
  }
  total_redun <- M1[1, l] ## total number of nodes in the redundant tree
  ## convert to a parent node vector
  node_redun <- rep(0, total_redun)
  for(i in 1:(total_redun-1)) {
    ind_ij <- which(M1 == i, arr.ind=TRUE)[1, ]
    ## the parent of node i is the node index to its right
    node_redun[i] <- M1[ind_ij[1], ind_ij[2]+1]
  }

  ## A trimmed/standardized version, parent node vector without redundancy
  M2 <- M ## node index matrix, index goes from 1 to # of nodes
  M2[, 1] <- 1:p
  for (j in 2:l) {
    M2[, j] <- match(M[, j], unique(M[, j]), nomatch = 0)
    tempseq <- M2[, j]
    for (k in max(tempseq):1) {
      if (length(unique(M2[which(tempseq == k), j-1])) == 1) {
        tempseq[tempseq == k] <- 0
        tempseq[tempseq > k] <- tempseq[tempseq > k] - 1
      }
    }
    tempseq[tempseq != 0] <- tempseq[tempseq != 0] + max(M2[, j-1])  ## unique new node index
    tempseq[tempseq == 0] <- M2[tempseq == 0, j-1]  ## carry over
    M2[, j] <- tempseq
  }

  total <- M2[1, l] ## total number of nodes in the redundant tree
  ## convert to a parent node vector
  node <- rep(0, total)
  for(i in 1:(total-1)) {
    ind_ij <- tail(which(M2 == i, arr.ind=TRUE), 1)
    ## the parent of node i is the node index to its right
    node[i] <- M2[ind_ij[1], ind_ij[2]+1]
  }
  return(list(node = node, M2 = M2))
}
