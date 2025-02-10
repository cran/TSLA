# require libraries
#library(ape)
#library(phytools)
#library(data.tree)

#' Generate aggregated features
#'
#' Function that generates aggregated features based on the TSLA output.
#'
#' @param TSLA.object A fit output from \code{TSLA.fit()},
#' or the \code{TSLA.fit} object in \code{cv.TSLA()}.
#' @param X_2 Expanded design matrix in matrix form.
#' @param X_2.org Original design matrix in matrix form.
#' @param lambda.index Index of the \eqn{\lambda} value selected.
#' @param alpha.index Index of the \eqn{\alpha} value selected. The \eqn{\alpha} is
#' the tuning parameter for generalized lasso penalty.
#' @return A data.frame of the aggregated feature.
#' \item{dataset}{aggregated features.}
#' @export
#'
getaggr <- function(TSLA.object, X_2, X_2.org,
                    lambda.index, alpha.index){
  root <- 1
  # gamma: final estimation of gamma coefficient
  # A: store A matrix from TSLA.object
  # ext_x X_2: extended x matrix file
  # splitfile: contain training obs id
  # g_idx: g_idx from TSLA.object
  # treeplot: M2 matrix
  # X_2.org: original design matrix
  # root: whether the root node is 0 or not; 0: zero, 1: not zero
  A <- TSLA.object$tree.object$A
  g_idx <- TSLA.object$tree.object$g_idx
  treenodes <- TSLA.object$tree.object$M2
  rmid <- TSLA.object$rmid
  gamma <- TSLA.object$gammacoef[ , lambda.index, alpha.index]
  groupnorm <- TSLA.object$groupnorm[ , lambda.index, alpha.index]

  rmid <- rmid[rmid <= ncol(X_2.org)]
  org_node <- treenodes[1:(ncol(X_2.org)-length(rmid)), ]
  org_node <- as.data.frame(org_node)
  if(length(rmid) >0 ){
    org_node$nodename <- colnames(X_2.org)[-rmid]
  }else{
    org_node$nodename <- colnames(X_2.org)
  }
  p <- nrow(treenodes) # number of leaves in expanded tree
  nall <- max(treenodes) # id of the root node
  groupsid <- list()
  for(i in 1:(nall-p)){
    # exclude nodename column
    idx <- which(org_node[, 1:(ncol(org_node)-1)] == i+p, arr.ind = TRUE)
    idx2 <- idx[, 2]-1
    # get direct child node id
    groupsid[[i]] <- setdiff(unique(unlist(org_node[idx[, 1], idx2])),i+p)
  }

  idg <- which(groupnorm/sum(groupnorm) > 1/nrow(g_idx))
  ##################################################################3
  ## get aggregation structure
  iduse <- setdiff(1:ncol(X_2.org), rmid)
  if(length(rmid) > 0){
    x <- X_2.org[, iduse]
  }else{
    x <- X_2.org
  }
  if(root == 1){
    ngroup <- nrow(g_idx)-1
    ninter <- ncol(A)
  }else{
    ngroup <- nrow(g_idx)
    ninter <- ncol(A)+1
  }
  internalid <- (nrow(A)+1):ninter
  xinter <- matrix(0, nrow = nrow(X_2.org), ncol = ninter)
  xinter <- data.frame(xinter)
  xinter[, 1:ncol(x)] <- x
  xaggr <- matrix(1, nrow = nrow(x), ncol=1)

  for(i in 1:ngroup){
    childnode <- groupsid[[i]]
    if(!(i %in% idg)){ # zero groups
      ifaggr <- apply(xinter[, childnode], 2, sum)
      nifaggr <- sum(ifaggr > 0)
      leadnodeid <- internalid[i]
      if( nifaggr > 0 ){
        if(nifaggr == 1){
          xinter[, leadnodeid] <- xinter[, childnode[ifaggr > 0 ]]
        }else{
          xinter[, leadnodeid] <-
            as.numeric(apply(xinter[, childnode[ifaggr > 0 ]], 1, sum)>0)
        }
      }
    }else{ # nonzero groups
      ifaggr <- apply(xinter[, childnode], 2, sum)
      #childnode <- childnode[ifaggr > 0 ]
      xaggr <- cbind(xaggr, xinter[, childnode[ifaggr > 0 ]])
    }

  }
  if(root == 1 & nrow(g_idx) %in% idg){
    if(sum(xinter[, ninter])>0){
      xaggr <- cbind(xaggr, xinter[, ninter])
    }
  }
  if(ncol(xaggr) == 1){
    if(root == 0){
      nofeature <- 1
      dataset <- NULL
      #warning('warning: no available structure!')
    }else{
      if(nrow(g_idx) %in% idg){
        dataset <- data.frame(x = as.numeric(apply(x, 1, sum)>0))
        nofeature <- 0
      }else{
        nofeature <- 1
        dataset <- NULL
        #warning('warning: no available structure!')
      }
    }
  }else{
    dataset <- data.frame(xaggr[, 2:ncol(xaggr)])
    nofeature <- 0
  }
  return(dataset)
}



#' Plot aggregated structure
#'
#' Return a tree plot.
#'
#' @param TSLA.object A fit output from \code{TSLA.fit()},
#' or the \code{TSLA.fit} object in \code{cv.TSLA()}.
#' @param X_2 Expanded design matrix in matrix form.
#' @param X_2.org Original design matrix in matrix form.
#' @param lambda.index Index of the \eqn{\lambda} value selected.
#' @param alpha.index Index of the \eqn{\alpha} value selected. The \eqn{\alpha} is
#' the tuning parameter for generalized lasso penalty.
#' @return A plot
#' @export
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
#' plot_TSLA(simu.cv$TSLA.fit, x2, x2.org, simu.cv$lambda.min.index, simu.cv$alpha.min.index)
#'
#'
#'
plot_TSLA <- function(TSLA.object, X_2, X_2.org,
                      lambda.index, alpha.index){
  root <- 1
  A <- TSLA.object$tree.object$A
  g_idx <- TSLA.object$tree.object$g_idx
  treenodes <- TSLA.object$tree.object$M2
  rmid <- TSLA.object$rmid
  gamma <- TSLA.object$gammacoef[ , lambda.index, alpha.index]
  groupnorm <- TSLA.object$groupnorm[ , lambda.index, alpha.index]
  rmid <- rmid[rmid <= ncol(X_2.org)]
  org_tree <- as.data.frame(treenodes[1:(ncol(X_2.org)-length(rmid)), ])
  org_node <- org_tree
  org_node <- as.data.frame(org_node)
  if(length(rmid) >0 ){
    org_node$nodename <- colnames(X_2.org)[-rmid]
    org_tree[, 1] <- colnames(X_2.org)[-rmid]
  }else{
    org_node$nodename <- colnames(X_2.org)
    org_tree[, 1] <- colnames(X_2.org)
  }

  # flip column order
  org_tree <- org_tree[, ncol(org_tree):1]
  # generate pathString for the as.Node function
  for(i in 1:(ncol(org_tree)-1)){
    types <- table(org_tree[, i])
    singleton <- which(types == 1)
    org_tree[which(org_tree[, i] %in% names(types)[singleton]), i:(ncol(org_tree)-1)] <- NA
    org_tree[which(org_tree[, i] == org_tree[, i+1]), i] <- NA
  }
  org_tree$pathString <- apply(org_tree[,1:ncol(org_tree)], 1,paste, collapse = '/')
  org_tree$pathString <- gsub(pattern = "/NA", "", org_tree$pathString)
  org_tree$pathString <- gsub(pattern = " ", "", org_tree$pathString)
  ############## full tree plot(no aggregation pattern )
  # full tree structure before any aggregation
  tree <- as.Node(org_tree)
  tree <- as.phylo(tree)
  tree$edge.length[1:length(tree$edge.length)] <- 100 / 6
  # plot(tree, direction = "downwards",
  #      adj = 0.5,
  #      lwd=0.4, label.offset = 1,
  #      font = 0.9,
  #      cex = 0.7, srt = -180, cex.main=0.9)
  # ## make leaf node to be gray
  # tiplabels(pch = 19, col = "grey", cex = 0.8,
  #           # bg = RColorBrewer::brewer.pal(n = 8, "Set1")[6])
  #           bg = "white")
  # ## make internal node to be grey
  # nodelabels(pch = 19, col = c("grey"), cex = 0.8,
  #            # bg = RColorBrewer::brewer.pal(n = 8, "Set1")[6])
  #            bg = "white")

  ################################################################
  ################################################################
  ## find lead node for each group
  p <- nrow(treenodes)
  nall <- max(treenodes)
  groups <- list()
  groupsid <- list()
  groupsidall <- list()
  for(i in 1:(nall-p)){
    idx <- which(org_node[, 1:(ncol(org_node)-1)] == i+p, arr.ind = TRUE)
    # descendant node name in each group: only for main effect nodes
    groups[[i]] <- unique(org_node$nodename[idx[,1]])
    # child node id
    groupsid[[i]] <- setdiff(unique(unlist(org_node[idx[, 1], idx[, 2]-1])),i+p)
    # child node id for both main and interaction nodes
    idxall <-  which(treenodes[, 1:(ncol(org_node)-1)] == i+p, arr.ind = TRUE)
    groupsidall[[i]] <- setdiff(unique(unlist(treenodes[idxall[, 1], idxall[, 2]-1])), i+p)
  }
  # get all descendant of a lead node
  leadnode <- (p+1):nall # lead node id
  groupsiddes <- groupsid
  for(i in 1:(nall-p)){
    contains <- which(leadnode %in% groupsiddes[[i]])
    if(sum(contains)>0){
      # case when some child nodes are internal nodes
      groupsiddes[[i]] <- unique(unlist(groupsiddes[contains]))
      groupsiddes[[i]] <- unique(c(unlist(groupsiddes[[i]]), groupsid[[i]]))
    }else{
      # case when all child nodes are leaf nodes
      groupsiddes[[i]] <- groupsiddes[[i]]
    }
  }
  # idg: nonzero group index in group order
  idg <- which(groupnorm / sum(groupnorm) > 1 / nrow(g_idx))
  lead_node_nonzero <- p+idg # this is nonzero groups's lead node id

  ###################################################
  ## discuss root node
  # if(nall-p == nrow(gidx)){
  #   root <- 0
  # }else{
  #   root <- 1
  # }

  ############################################################
  ##############################################################
  ## make graph
  ## let some node to be white
  ## make all leaf node to be open: pch = 21
  tips.pch <- rep(21, nrow(org_node))
  ## find all leaf nodes of nonzero groups
  node_nonzero <- NULL
  all_nonzero_node <- unique(unlist(groupsid[idg]))
  all_nonzero_leafnode <- intersect(all_nonzero_node, 1:nrow(org_node))
  all_nonzero_leafnode_name <- org_node$nodename[all_nonzero_leafnode]
  # get direct parent of all leaf nodes
  #leafnode_parent <-
  #  apply(org_node[,1:(ncol(org_node)-1)], 1,
  #        function(x) min(x[x != x[1]]))
  #leaf_node_nonzero <- which(leafnode_parent %in% lead_node_nonzero)
  #leaf_node_nonzero_name <- org_node$nodename[leaf_node_nonzero]
  #let nonzero leaf node to be solid
  tips.pch[which(tree$tip.label %in% all_nonzero_leafnode_name)] <- 19


  if(root == 1){
    if(nall %in% (idg+p-1)){
      sparse_leaf_node <- NULL
    }else{
      sparse_leaf_node <- apply(org_node[, 1:(ncol(org_node)-1)], 1,
                                FUN=function(x){sum(x %in% lead_node_nonzero)})
    }
  }else{
    sparse_leaf_node <- apply(org_node[, 1:(ncol(org_node)-1)], 1,
                              FUN=function(x){sum(x %in% lead_node_nonzero)})
  }
  # all ancester and itself are zero nodes
  sparse_leaf_node_name <- org_node$nodename[which(sparse_leaf_node == 0)]
  tips.pch[which(tree$tip.label %in% sparse_leaf_node_name)] <- 4 # cross sign

  #############################################################3
  ## find internal nodes of nonzero groups
  inter_node_nonzero <- unique(unlist(groupsid[idg]))
  inter_node_nonzero <- setdiff(inter_node_nonzero, 1:nrow(treenodes))

  node.pch <- rep(21, length(tree$node.label))
  node.pch[which(tree$node.label %in% inter_node_nonzero)] <-19

  ## find internal nodes with nonzero descendant
  des_nonzero <- lapply(groupsiddes, FUN=function(x){sum(x %in% c(all_nonzero_leafnode, inter_node_nonzero))})
  des_nonzero <- unlist(des_nonzero)
  inter_node_nonzero_des <- which(des_nonzero > 0) + p
  node.pch[which(tree$node.label %in% inter_node_nonzero_des)] <-19
  ## discuss root node
  if(root == 1){
    if(nall %in% (p+idg-1)){
      node.pch[which(tree$node.label == nall)] <- 19
    }
  }

  ###########################################################3
  ## check which edge to be dashed: 1 is solid, 2 is dashed
  ## tree edge index has different ordering system
  # relate these two system
  idtrans <- data.frame(old = (nrow(treenodes)+1):nall, new = nrow(org_node)+1)
  id <- nrow(org_node) + 1
  idset <- c(nall)
  for(i in 1:nrow(org_node)){
    # start from V5
    for(j in (ncol(org_node)-2):2){
      idold <- org_node[i, j]
      if((!(idold %in% idset)) & idold > nrow(org_node)){
        id <- id + 1
        idset <- c(idset, idold)
        idtrans$new[idtrans$old == idold] <- id
      }
    }
  }

  edge.type <- rep(1, length(tree$edge.length))
  # dashed line when point to zero leaf nodes
  edge.type[which(apply(tree$edge, 1,
                        function(x) {
                          as.numeric(x[2] %in% which(tips.pch %in% c(4, 21)))
                        })==1)] <- 2
  # dashed line when point to zero internal nodes (exclude lead node of selected groups)
  zero_nlabel <- idtrans$new[!idtrans$old %in% inter_node_nonzero]
  select_lead_node <- idtrans$new[idtrans$old %in% (idg+nrow(treenodes))]
  zero_nlabel <- setdiff(zero_nlabel, select_lead_node)
  edge.type[which(apply(tree$edge, 1,
                        function(x) {
                          as.numeric(x[2] %in% zero_nlabel)
                        })==1)] <- 2
  plot(tree, direction = "rightward",
       adj = 0.1,
       lwd=0.4, label.offset = 1,
       font = 0.5,
       cex = 0.7, srt = 0,
       edge.width = 3 - edge.type,
       edge.lty = edge.type)


  tiplabels(pch = tips.pch, col = "black", cex = 0.8,
            bg = "white")
  nodelabels(pch = node.pch, col = c("black"), cex = 0.8,
            bg = "white")

}
