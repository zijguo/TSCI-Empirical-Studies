library(ranger)
library(Matrix)


# Random Forest with data splitting
#
# param D continuous or binary, n by 1 treatment vector.
# param Z continuous or binary, n by 1 Instrumental Variable, only one instrument is implemented for violation space selection now.
# param X continuous or binary, n by p_x baseline covariates matrix.
# param num.trees number of trees.
# param mtry number of covariates to possibly split at in each node.
# param max.depth maximal depth of each tree.
# param min.node.size minimal size of each leaf node.
# param split.prop a value between 0 and 1, the proportion of samples we use in A1.
# param MSE.thol a very large value of MSE, default by 1e12, used for the start of hyper-parameter selection.
# param forest.save save the Random Forest output or not, default by TRUE.
#
# @return:
#     \item{\code{forest.A2}}{random forest built on subsample A2, available if forest.save=TRUE.}
#     \item{\code{params.A2}}{best hyper-parameters for forest.A2 selected by out-of-bag error.}
#     \item{\code{A1.ind}}{indices of subsample A1.}
#     \item{\code{nodes.A1}}{a n_A1 by num.trees matrix, rows refer to different samples, columns refer to different trees, the entrees are leaf node indices of each sample in each tree.}
#     \item{\code{MSE.oob}}{minimal out-of-bag error using the best hyper-parameters.}
#
#
TSCI.RF.fit <- function(D,Z,X,num.trees,mtry,max.depth,min.node.size,split.prop,MSE.thol=1e12,forest.save=FALSE) {
  W <- as.matrix(cbind(Z,X)); D <- as.matrix(D)
  n <- NROW(W); p <- NCOL(W)
  Data <- data.frame(cbind(D,W))
  names(Data) <- c("D", paste("W", 1:p, sep = ""))
  # grid search
  params.grid <- expand.grid(
    num.trees = num.trees,
    mtry = mtry,
    max.depth = max.depth,
    min.node.size = min.node.size
  )
  # split the data into two parts A1 and A2
  # use A2 to build the random forest and use A1 to predict
  n.A1 <- round(split.prop*n)
  A1.ind <- 1:n.A1
  Data.A1 <- Data[A1.ind,]
  Data.A2 <- Data[-A1.ind,]
  forest.A2 <- NULL;
  MSE.oob.A2 <- MSE.thol
  params.A2 <- NULL
  ### use oob error to do hyper-parameter tuning
  for (i in 1:nrow(params.grid)) {
    temp.A2 <- ranger(D~., data = Data.A2,
                      num.trees=params.grid$num.trees[i],
                      mtry=params.grid$mtry[i],
                      max.depth = params.grid$max.depth[i],
                      min.node.size = params.grid$min.node.size[i]
    )
    if (temp.A2$prediction.error <= MSE.oob.A2) {
      forest.A2 <- temp.A2
      params.A2 <- params.grid[i,]
      MSE.oob.A2 <- temp.A2$prediction.error
    }
  }
  # leaf nodes information of A1 on the Random Forest built on A2
  nodes.A1 <- predict(forest.A2, data = Data.A1, type = "terminalNodes")$predictions
  returnList <- list(forest.A2 = forest.A2,
                     params.A2 = params.A2,
                     A1.ind = A1.ind,
                     nodes.A1 = nodes.A1,
                     MSE.oob.A2 = MSE.oob.A2)
  if (!forest.save) returnList <- returnList[-1]
  returnList
}


# Weight matrix computation of random forest
#
# param: nodes A n_A1 by num.trees matrix, rows refer to different samples, columns refer to different trees, the entrees are leaf node indices of each sample in each tree.
#
# return: \item{\code{out.weight}}{A n_A1 by n_A1 symmetric sparse weight matrix of class dgCMatrix with the ith row represents the weights of outcome of each sample on the prediction of the ith outcome.}
TSCI.RF.weight <- function(nodes) {
  n.A1 <- NROW(nodes); num.trees <- NCOL(nodes)
  out.weight <- matrix(0,n.A1,n.A1)
  for (j in 1:num.trees) {
    weight.mat <- matrix(0,n.A1,n.A1) # weight matrix for single tree
    unique.nodes <- unique(nodes[,j])
    for (i in 1:length(unique.nodes)) {
      ind <- nodes[,j]==unique.nodes[i] # indices of samples in the node
      num.samples <- sum(ind) # number of samples in the node
      w <- 1/(num.samples-1)  # weight, to remove self-prediction
      weight.vec <- ifelse(ind,yes=w,no=0)
      weight.mat[ind,] <- matrix(rep(weight.vec,num.samples),num.samples,byrow=T)/num.trees
    }
    diag(weight.mat) <- 0 # remove self prediction
    out.weight <- out.weight + weight.mat
  }
  out.weight <- Matrix(out.weight, sparse = T) # sparse matrix to save memory
  return(out.weight)
}








### naiveRF.stat
### Function: Compute the plug-in TSLS using random forest with
###           data splitting as the first stage
### Input: Y: continuous, n by 1 outcome vector
###        D: continuous or binary, n by 1 treatment vector
###        VW: the instruments-covariates matrix
###        A1.ind: a set of indices indicating the samples in A1
###        weight: a n by n symmetric sparse weight matrix(class dgCMatrix), the ith row
###                represents the weights of each sample on the prediction of the ith outcome
### Output: betaHat: the estimated treatment effect
###         sd: the standard deviation of betaHat
###         singularity: logical, whether t(VW.rep)%*%VW.rep is singular
naiveRF.stat <- function(Y, D, VW, A1.ind, weight, lam=0.05) {
  n.A1 <- length(A1.ind); r.VW <- ncol(VW)
  Y.A1 <- Y[A1.ind]; VW.A1 <- VW[A1.ind,]; D.A1 <- D[A1.ind]
  D.rep <- as.matrix(weight%*%D.A1)
  
  ### point estimator
  reg.rf <- lm(Y.A1~D.rep+VW.A1)
  betaHat <- coef(reg.rf)[2]
  
  ### noise level of outcome model
  SigmaSqY <- get.sigma(betaHat, Y, D, VW)
  
  D.resid <- resid(lm(D.rep~VW.A1))
  sd <- sqrt(SigmaSqY/sum(D.resid^2))
  
  returnList <- list(betaHat = betaHat,
                     sd = sd,
                     SigmaSqY = SigmaSqY)
  returnList
}


### TSRF.full
### Function: Model fitting function for TWo Stage Random Forest using
###           full data. Use out-of-bag error for hyper-parameter tuning
### Input: X: continuous or binary, n by p_x covariates
###        Y: continuous or binary, n by 1 outcome vector
###        k: integer, number of subsamples for cross fitting
###        num.trees: integer, the number of trees in random forest
###        mtry: integer, the number of covariates to split at each node
###        max.depth: integer, the maximal depth of each tree, 0 refers to unlimited depth
###        min.node.size: integer, the minimal size(# samples in it) of each leaf node
###        MSE.thol: numeric, a large value of MSE, used for the start of hyper-parameter selection
###        forest.save: logic, to save the random forest object or not, default by FALSE to save memory
### Output: forest: random forest object using full data, available if forest.save=TRUE
###         params: a list of best hyper-parameters selected by the out-of-bag error
###         predicted.values: the predicted values of outcome Y using full data
###         nodes: a n by num.trees nodes information matrix, similar to nodes.A1 in TSRF.fit
###         MSE.oob: the minimal out-of-bag error using the best hyper-parameters
TSRF.full <- function(X,Y,num.trees=200,mtry=NULL,max.depth=0,min.node.size=5,MSE.thol=1e6,forest.save=F) {
  n <- nrow(X);p<-ncol(X)
  if (is.null(mtry)) mtry <- round(p/3)
  
  Data <- data.frame(cbind(Y, X))
  names(Data) <- c("Y", paste("X", 1:p, sep = ""))
  ### grid search
  params.grid <- expand.grid(
    num.trees = num.trees,
    mtry = mtry,
    max.depth = max.depth,
    min.node.size = min.node.size
  )
  
  # use oob error to do hyper-parameter tuning
  forest <- params <- NA
  MSE.oob <- MSE.thol
  for (i in 1:nrow(params.grid)) {
    temp<- ranger(Y~., data = Data,
                  num.trees=params.grid$num.trees[i],
                  mtry=params.grid$mtry[i],
                  max.depth = params.grid$max.depth[i],
                  min.node.size = params.grid$min.node.size[i]
    )
    if (temp$prediction.error < MSE.oob) {
      forest <- temp
      params <- params.grid[i,]
      MSE.oob <- temp$prediction.error
    }
  }
  
  # conduct prediction
  predicted.values <- predict(forest,data = Data, type="response")$predictions
  nodes <- predict(forest,data = Data, type="terminalNodes")$predictions
  
  returnList = list(forest = forest,
                    params = params,
                    predicted.values = predicted.values,
                    nodes = nodes,
                    MSE.oob = MSE.oob)
  if (!forest.save) returnList <- returnList[-1]
  returnList
}


### TSRF.stat.full
### Function: Compute the Two Stage Random Forest using full data
### Input: Y: continuous, n by 1 outcome vector
###        D: continuous or binary, n by 1 treatment vector
###        VW: the instruments-covariates matrix
###        weight: a n by n symmetric sparse weight matrix(class dgCMatrix), the ith row
###                represents the weights of each sample on the prediction of the ith outcome
### Output: betaHat: the estimated treatment effect
###         sd: the standard deviation of betaHat
###         singularity: logical, whether t(VW.rep)%*%VW.rep is singular
TSRF.stat.full <- function(Y, D, VW, weight) {
  n <- length(Y); r.VW <- ncol(VW)
  D.rep <- as.matrix(weight%*%D); Y.rep <- as.matrix(weight%*%Y)
  VW.rep <- as.matrix(weight%*%VW)
  SigmaSqD = mean((D.rep-D)^2)
  
  ### point estimator
  reg.rf <- lm(Y.rep~D.rep+VW.rep)
  betaHat <- coef(reg.rf)[2]
  
  ### noise level of outcome model
  SigmaSqY <- get.sigma(betaHat, Y, D, VW)
  
  ### the standard error of the estimator
  D.resid <- resid(lm(D.rep~VW.rep))
  iv.str <- sum(D.resid^2)
  sd <- sqrt(SigmaSqY*(t(D.resid)%*%as.matrix(weight)%*%as.matrix(weight)%*%D.resid)/(iv.str^2))
  
  returnList <- list(betaHat = betaHat,
                     sd = sd,
                     SigmaSqY = SigmaSqY,
                     iv.str = iv.str,
                     SigmaSqD = SigmaSqD)
  returnList
}



### get.sigma
### Function: a helper function for the estimate of noise level
### Input: betaHat: the estimated treatment effect by data splitting estimator
###        Y: continuous, n by 1 outcome vector
###        D: continuous or binary, n by 1 treatment vector
###        VW: the instruments-covariance matrix expanded by corresponding basis,
###            first column is intercept
### Output: the estimate of noise level, typically SigmaSqY
get.sigma <- function(betaHat, Y, D, VW) {
  n <- length(Y)
  Y.D <- Y-betaHat*D
  Y.resid <- resid(lm(Y.D~VW))
  sum(Y.resid^2)/(n-NCOL(VW)-1)
}



