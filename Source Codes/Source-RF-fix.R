library(Matrix)
library(ranger)


#' Random Forest with data splitting
#'
#' @param D continuous or binary, n by 1 treatment vector.
#' @param Z continuous or binary, n by 1 Instrumental Variable, only one instrument is implemented for violation space selection now.
#' @param X continuous or binary, n by p_x baseline covariates matrix.
#' @param num.trees number of trees.
#' @param mtry number of covariates to possibly split at in each node.
#' @param max.depth maximal depth of each tree.
#' @param min.node.size minimal size of each leaf node.
#' @param split.prop a value between 0 and 1, the proportion of samples we use in A1.
#' @param MSE.thol a very large value of MSE, default by 1e12, used for the start of hyper-parameter selection.
#' @param forest.save save the Random Forest output or not, default by TRUE.
#'
#' @return
#'     \item{\code{forest.A2}}{random forest built on subsample A2, available if forest.save=TRUE.}
#'     \item{\code{params.A2}}{best hyper-parameters for forest.A2 selected by out-of-bag error.}
#'     \item{\code{A1.ind}}{indices of subsample A1.}
#'     \item{\code{nodes.A1}}{a n_A1 by num.trees matrix, rows refer to different samples, columns refer to different trees, the entrees are leaf node indices of each sample in each tree.}
#'     \item{\code{MSE.oob}}{minimal out-of-bag error using the best hyper-parameters.}
#' @export
#'
#' @examples
#' \donttest{
#'
#' }
#'
#'
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


#' Weight matrix computation of random forest
#'
#' @param nodes A n_A1 by num.trees matrix, rows refer to different samples, columns refer to different trees, the entrees are leaf node indices of each sample in each tree.
#'
#' @return
#'     \item{\code{out.weight}}{A n_A1 by n_A1 symmetric sparse weight matrix of class dgCMatrix with the ith row represents the weights of outcome of each sample on the prediction of the ith outcome.}
#' @export
#'
#' @examples
#' \donttest{
#'
#' }
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




### TSCI.RF.fix
### Function: Violation space selection of Two Stage Curvature Identification using Random Forest
### Input: Y: continuous, n by 1 outcome vector
###        D: continuous or binary, n by 1 treatment vector
###        Cov.aug: Augmented combination of violation space and baseline covariates
###        A1.ind: Indices of samples in A1
###        weight: n.A1 by n.A1 weight matrix
###        Q: Number of violation space, including no violation 
###        intercept: Include intercept in the outcome model or not
### Output: produce the same outcome of function TSCI.RF
TSCI.RF.fix <- function(Y, D, Cov.aug, A1.ind, weight, Q, intercept=TRUE) {
  Y <- as.matrix(Y); D <- as.matrix(D); Cov.aug <- as.matrix(Cov.aug)
  # constants
  n <- length(Y); n.A1 <- length(A1.ind);
  Y.A1 <- Y[A1.ind]; D.A1 <- D[A1.ind]; Cov.aug.A1 <- Cov.aug[A1.ind,]
  # compute the representations, fix Y, V, W
  D.rep <- as.matrix(weight%*%D.A1)
  # the noise level of treatment model
  SigmaSqD <- mean((D.rep-D.A1)^2)
  # save estimates for selection part
  names <- paste("RF-fix-q",0:(Q-1),sep="")
  Coef.vec <- sd.vec <- rep(NA,Q)
  names(Coef.vec) <- names(sd.vec) <- names
  # the noise level of outcome model and covariance of epsilon and delta
  SigmaSqY <- SigmaYD <- denom <- numer <-  rep(NA,Q)
  names(SigmaSqY) <- names(SigmaYD) <- names(denom) <- names(numer) <- names
  
  ### fixed violation space, compute necessary inputs of selection part
  for (index in 1:Q) {
    q <- index-1
    if (q==Q-1) {
      if (intercept) {
        Cov.aug.A1.q = cbind(1,Cov.aug.A1)
        Cov.aug.q = cbind(1,Cov.aug)
      }
    } else {
      Cov.aug.q = Cov.aug[,-(1:(Q-1-q))]
      Cov.aug.A1.q = Cov.aug.A1[,-(1:(Q-1-q))]                   
      if (intercept) {
        Cov.aug.A1.q = cbind(1,Cov.aug.A1.q)
        Cov.aug.q = cbind(1,Cov.aug.q)
      }
    }
    P = diag(1,n.A1,n.A1) - Cov.aug.A1.q%*%solve(t(Cov.aug.A1.q)%*%Cov.aug.A1.q)%*%t(Cov.aug.A1.q)
    D.RSS = as.numeric(t(D.A1)%*%P%*%D.rep) # denominator
    Cov.YD = as.numeric(t(Y.A1)%*%P%*%D.rep) # numerator
    trace.T = sum(diag(P%*%weight))
    SigmaYD[index] = sum((D.A1-D.rep)*resid(lm(Y.A1-D.A1*(Cov.YD/D.RSS)~Cov.aug.A1.q)))/(n.A1-ncol(Cov.aug.A1.q))
    Coef.vec[index] = (Cov.YD-SigmaYD[index]*trace.T)/D.RSS
    SigmaSqY[index] = get.sigma(Coef.vec[index],Y,D,Cov.aug.q)
    sd.vec[index] = sqrt(SigmaSqY[index]*t(D.rep)%*%P%*%D.rep)/D.RSS
    denom[index] = D.RSS
    numer[index] = Cov.YD
  }
  
  returnList = list(Coef.vec = Coef.vec,
                    sd.vec = sd.vec,
                    SigmaSqY = SigmaSqY,
                    SigmaSqD = SigmaSqD,
                    SigmaYD = SigmaYD,
                    denom = denom,
                    numer = numer
                    )
  returnList
}


### get.sigma
### Function: A helper function for the estimate of noise level
### Input: betaHat: the estimated treatment effect by data splitting estimator
###        Y: continuous, n by 1 outcome vector
###        D: continuous or binary, n by 1 treatment vector
###        Cov.aug: the instruments-covariance matrix expanded by corresponding basis,
###            first column is intercept
### Output: the estimate of noise level, typically SigmaSqY
get.sigma <- function(betaHat, Y, D, Cov.aug) {
  n <- length(Y)
  Y.D <- Y-betaHat*D
  Y.resid <- resid(lm(Y.D~Cov.aug))
  sum(Y.resid^2)/(n-NCOL(Cov.aug)-1)
}

