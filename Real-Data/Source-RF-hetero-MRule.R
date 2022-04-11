# This is for maimonide's rule and the difference is that we add violation forms to first stage's random forest fitting

library(Matrix)
library(ranger)


#' @title Two Stage Curvature Identification with Random Forest
#' @description Implementation of Two Stage Curvature Identification with random forest. An IV strength test and the violation space selection are implemented. It estimates the causal effect in the instrumental variable framework and construct confidence intervals
#'
#' @param Y continuous, n by 1 outcome vector
#' @param D continuous or binary, n by 1 treatment vector
#' @param Z continuous or binary, n by 1 instrumental variable vector
#' @param X continuous or binary, n by p baseline covariates matrix
#' @param intercept logic, including intercept or not in the outcome model, default by TRUE
#' @param vio.space a matrix with each column corresponding to a violation form of Z, default by polynomial violation space to the degree of 3
#' @param layer logic, layer selection of violation space or not, default by TRUE
#' @param split.prop a value between 0 and 1, the proportion of samples used for inference in data splitting. The rest of samples are used for building random forest
#' @param Omega a user-provided weight/transformation matrix, random forest will be performed if not provided. This must be a n.A1 by n.A1 matrix where n.A1 = round(n*split.prop), denoting sample size for inference
#' @param num.trees number of trees, default by 200
#' @param mtry number of covariates to possibly split at in each node of the tree, default by sequence from (p+1)/3 to 2(p+1)/3 for cross validation in hyper-parameter tuning
#' @param max.depth maximal depth of each tree, default by 0 referring to unlimited depth
#' @param min.node.size minimal size of each leaf node, default by the set of 5, 10, 15 for cross validation in hyper-parameter tuning
#' @param str.thol minimal value of the threshold of IV strength test, default by 20
#' @param alpha confidence level used to construct confidence intervals, default by 0.05
#'
#' @return
#'     \item{\code{Coef.vec}}{original and bias-corrected point estimators in fixed violation space}
#'     \item{\code{sd.vec}}{standard error of original and bias-corrected point estimators in fixed violation space}
#'     \item{\code{Coef.robust}}{original and bias-corrected point estimators in the violation space selected by the comparison and robust methods}
#'     \item{\code{sd.robust}}{standard errors of original and bias-corrected point estimators in the violation space selected by the comparison and robust methods}
#'     \item{\code{CI.robust}}{confidence intervals of original and bias-corrected estimators in the violation space selected by the comparison and robust methods}
#'     \item{\code{SigmaSqY}}{estimated variance of the error term in outcome model}
#'     \item{\code{SigmaSqD}}{estimated variance of the error term in the treatment model}
#'     \item{\code{SigmaSqY.Qmax}}{esimated variance of the error term in the outcome model, in the largest violation space that the IV strength test is passed}
#'     \item{\code{iv.str}}{IV strength in fixed violation space}
#'     \item{\code{iv.thol}}{IV strength test threshold in fixed violation space}
#'     \item{\code{Qmax}}{index of largest violation space selected by IV strength test}
#'     \item{\code{q.comp}}{index of estimated violation space using comparison method}
#'     \item{\code{q.robust}}{index of estimated violation space using robust method}
#'     \item{\code{invalidity}}{invalidity of TSLS, at least linear violation or violation in the first violation form if TRUE}
#' @export
#'
#' @examples
#' \dontrun{
#' ### dimension
#' p = 10
#' ###
#' n = 100
#' ### interaction value
#' inter.val = 1
#' ### a denotes the IV strength
#' a = 1
#' ### violation strength
#' tau = 1
#' f <- function(x){a*(1*sin(2*pi*x) + 1.5*cos(2*pi*x))}
#' rho1=0.5
#' Cov<-(A1gen(rho1,p+1))
#' mu<-rep(0,p+1)
#' beta=1
#' alpha<-as.matrix(rep(-0.3,p))
#' gamma<-as.matrix(rep(0.2,p))
#' inter<-as.matrix(c(rep(inter.val,5),rep(0,p-5)))
#'
#'
#' ### generate the data
#' mu.error<-rep(0,2)
#' Cov.error<-matrix(c(1,0.5,0.5,1),2,2)
#' Error<-mvrnorm(n, mu.error, Cov.error)
#' W.original<-mvrnorm(n, mu, Cov)
#' W<-pnorm(W.original)
#' Z<-W[,1]
#' X<-W[,-1]
#' ### generate the treatment variable D
#' D=f(Z)+X%*%alpha+Z*X%*%inter+Error[,1]
#' ### generate the outcome variable Y
#' Y=D*beta+ tau*Z+ X%*% gamma+Error[,2]
#'
#'
#' ### random forest
#' output.RF <- TSCI.RF(Y,D,Z,X)
#' # point estimate
#' output.RF$Coef.robust
#' # standard error
#' output.RF$sd.robust
#' # confidence interval
#' output.RF$CI.robust
#' }
#'
#' @importFrom MASS mvrnorm
#' @importFrom Matrix Matrix
#' @importFrom stats coef lm poly predict qnorm quantile resid rnorm
#' @import ranger
#'
TSCI.RF <- function(Y,D,Z,X,intercept=TRUE,vio.space=NULL,layer=TRUE,split.prop=2/3, Omega = NULL,
                    num.trees=NULL,mtry=NULL,max.depth=NULL,min.node.size=NULL,str.thol=10,alpha=0.05) {
  Y <- as.matrix(Y); D <- as.matrix(D); Z <- as.matrix(Z); X <- as.matrix(X);
  # constants
  n <- NROW(X); p <- NCOL(X) + NCOL(Z)
  # default value for hyper-parameters
  if (is.null(num.trees)) num.trees <- 200
  if (is.null(mtry)) mtry <- seq(round(p/3),round(2*p/3),by=1)
  if (is.null(max.depth)) max.depth <- 0
  if(is.null(min.node.size)) min.node.size <- c(5,10,20)
  # define the vio.space if not specified
  if (is.null(vio.space)) {
    Q = 4
    vio.space <- matrix(NA,nrow(Z),0)
    for (q in 1:(Q-1)) {
      vio.space <- cbind(Z^q,vio.space)
    }
  } else {
    Q = NCOL(vio.space) + 1
  }
  # define the augmentation of covariates,
  # which is the combination of violation space and baseline covariates
  Cov.aug <- cbind(vio.space,X)
  n.A1 <- round(split.prop*n)
  A1.ind <- 1:n.A1
  
  if (is.null(Omega)) {
    # Treatment model fitting
    forest <- TSCI.RF.fit(D,Cov.aug,num.trees=num.trees,mtry=mtry,max.depth=max.depth,min.node.size=min.node.size,split.prop=split.prop)
    # Compute weight matrix
    weight <- TSCI.RF.weight(forest$nodes.A1)
  } else {
    weight = Omega
  }
  
  # Selection
  outputs <- TSCI.RF.Selection(Y,D,Cov.aug,A1.ind,weight=weight,Q=Q,intercept=intercept,layer=layer,str.thol=str.thol,alpha=alpha)
  return(outputs)
}


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
TSCI.RF.fit <- function(D,W,num.trees,mtry,max.depth,min.node.size,split.prop,MSE.thol=1e12,forest.save=FALSE) {
  W <- as.matrix(W); D <- as.matrix(D)
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


###### helper functions for violation space selection

### TSCI.RF.stat
### Function: Compute the necessary statistics for TSCI with Random Forest
### Input: Y.rep: continuous, a n.A1 by 1 vector denoting the outcome representation in subsample A1
###        D.rep: continuous or binary, a n.A1 by 1 vector denoting the treatment representation in subsample A1
###        Cov.rep: continuous or binary, the augmented instruments-covariates matrix representation in subsample A1
###        betaHat: Estimated treatment effect by data splitting estimator
###        weight: n.A1 by n.A1 weight matrix
###        n: full sample size
###        SigmaSqY: Estimated noise level of outcome model
###        SigmaSqD: Estimated noise level of treatment model
###        SigmaYD: Estimated covariance between error terms in outcome and treatment models
###        str.thol: the minimal value of the threshold of IV strength test, default by 20
### Output: betaHat: Estimated treatment effect by data splitting estimator
###         sd: the estimated standard deviation of data splitting estimator
###         betaHat.cor: the bias corrected estimator
###         sd.cor: the estimated standard deviation of bias corrected estimator
###         D.resid: Residuals of second stage regression
###         iv.str: IV strength
###         iv.thol: IV Strength Test threshold
###         explained.iv: t(D.A1)%*%T.V%*%T.V%*%D.A1, saved for violation space selection
TSCI.RF.stat <- function(D.rep, Cov.rep, weight, n, eps.hat, delta.hat, str.thol) {
  n.A1 <- length(D.rep); r.aug <- NCOL(Cov.rep)
  # compute the trace of T(V)
  # the trace of T matrix can be computed as RSS of each column of Omega on Cov.rep
  SigmaSqD = mean(delta.hat^2)
  RSS.V = rep(NA,n.A1)
  for (j in 1:n.A1) {
    RSS.V[j] <- sum(resid(lm(weight[,j]~Cov.rep))^2)
  }
  trace.T = sum(RSS.V)
  D.rep2 <- weight%*%D.rep
  D.resid <- resid(lm(D.rep~Cov.rep))
  D.RSS <- sum(D.resid^2)
  iv.str <- D.RSS/SigmaSqD
  # this is the numerator of the variance of betaHat
  explained.iv <- as.numeric(t(D.resid)%*%weight%*%weight%*%D.resid)
  sd <- sqrt(sum(eps.hat^2*(weight%*%D.resid)^2))/D.RSS
  
  
  
  ### standard errors of bias-corrected estimator
  # betaHat.cor <- betaHat - SigmaYD*trace.T/D.RSS
  # sd.cor <- sqrt((SigmaSqY*explained.iv+(SigmaSqD*SigmaSqY+SigmaYD^2)*(trace.T^2)/(n.A1-r.aug-1))/(D.RSS^2))
  
  # bootstrap for the threshold of IV strength test
  boot.vec <- rep(NA,300)
  delta.cent = delta.hat - mean(delta.hat)
  for (i in 1:300) {
    delta = rep(NA,n.A1)
    for (j in 1:n.A1) {
      U.j = rnorm(1)
      delta[j] = delta.cent[j]*U.j
    }
    
    delta.rep <- weight%*%delta
    delta.resid <- resid(lm(as.matrix(delta.rep)~Cov.rep))
    boot.vec[i] <- sum(delta.resid^2) + 2*sum(D.rep2*delta.resid)
  }
  iv.thol <- quantile(boot.vec,0.975)/SigmaSqD + max(2*trace.T, str.thol)
  # scale <- 1
  returnList <- list(
    # betaHat = betaHat,
    sd = sd,
    # betaHat.cor = betaHat.cor,
    # sd.cor = scale*sd.cor,
    D.resid = D.resid,
    iv.str = iv.str,
    iv.thol = iv.thol,
    explained.iv = explained.iv,
    trace.T = trace.T,
    RSS.V = RSS.V)
  returnList
}


### TSCI.RF.Selection
### Function: Violation space selection of Two Stage Curvature Identification using Random Forest
### Input: Y: continuous, n by 1 outcome vector
###        D: continuous or binary, n by 1 treatment vector
###        Cov.aug: Augmented combination of violation space and baseline covariates
###        A1.ind: Indices of samples in A1
###        weight: n.A1 by n.A1 weight matrix
###        Q: Number of violation space, including no violation
###        intercept: Include intercept in the outcome model or not
###        layer: logic, layer selection or not
###        str.thol: the minimal value of the threshold of IV strength test
### Output: produce the same outcome of function TSCI.RF
TSCI.RF.Selection <- function(Y, D, Cov.aug, A1.ind, weight, Q, intercept, layer, str.thol, alpha) {
  Y <- as.matrix(Y); D <- as.matrix(D); Cov.aug <- as.matrix(Cov.aug)
  # constants
  n <- length(Y); n.A1 <- length(A1.ind); r.aug <- NCOL(Cov.aug)
  Y.A1 <- Y[A1.ind]; D.A1 <- D[A1.ind]; Cov.aug.A1 <- Cov.aug[A1.ind,]
  # compute the representations
  Y.rep <- as.matrix(weight%*%Y.A1); D.rep <- as.matrix(weight%*%D.A1)
  Cov.rep <- as.matrix(weight%*%Cov.aug.A1)
  # the noise of treatment model
  delta.hat = D.A1 - D.rep
  SigmaSqD = mean(delta.hat^2)
  # save estimates for selection part
  names <- c(paste("RF-q",0:(Q-1),sep=""),paste("RF-Cor-q",0:(Q-1),sep=""))
  Coef.vec <- sd.vec <- rep(NA,2*Q)
  names(Coef.vec) <- names(sd.vec) <- names
  # IV strength test and signal strength test
  iv.str <- iv.thol <- rep(NA,Q)
  names(iv.str) <- names(iv.thol) <- paste("q",0:(Q-1),sep="")
  # the noise of outcome model
  eps.hat <- rep(list(NA),Q)
  # the numerator of variance
  explained.iv <- rep(NA,Q)
  names(explained.iv) <- paste("q",0:(Q-1),sep="")
  trace.T = explained.iv
  SigmaSqY = trace.T
  
  
  ### fixed violation space, compute necessary inputs of selection part
  # save D.resid for the computation of H and z.alpha
  D.resid <- RSS.vec <- rep(list(NA),Q)
  for (index in 1:Q) {
    q <- index-1
    if (q==Q-1) {
      if (intercept) {
        reg.rf <- lm(Y.rep~D.rep+Cov.rep)
        betaHat <- coef(reg.rf)[2]
      } else {
        reg.rf <- lm(Y.rep~D.rep+Cov.rep-1)
        betaHat <- coef(reg.rf)[1]
      }
      Coef.vec[index] <- betaHat
      eps.hat[[index]] = resid(lm(Y.A1-D.A1*betaHat~Cov.aug.A1))
      SigmaSqY[index] = mean(eps.hat[[index]]^2)
      stat.outputs <- TSCI.RF.stat(D.rep,Cov.rep,weight,n,eps.hat[[index]],delta.hat,str.thol=str.thol)
    } else {
      if (intercept) {
        reg.rf <- lm(Y.rep~D.rep+Cov.rep[,-(1:(Q-1-q))])
        betaHat <- coef(reg.rf)[2]
      } else {
        reg.rf <- lm(Y.rep~D.rep+Cov.rep[,-(1:(Q-1-q))]-1)
        betaHat <- coef(reg.rf)[1]
      }
      Coef.vec[index] <- betaHat
      eps.hat[[index]] = resid(lm(Y.A1-D.A1*betaHat~Cov.aug.A1[,-(1:(Q-1-q))]))
      SigmaSqY[index] = mean(eps.hat[[index]]^2)
      stat.outputs <- TSCI.RF.stat(D.rep,Cov.rep[,-(1:(Q-1-q))],weight,n,eps.hat[[index]],delta.hat,str.thol=str.thol)
    }
    # the statistics
    # Coef.vec[index+Q] <- stat.outputs$betaHat.cor
    sd.vec[index] <- stat.outputs$sd
    sd.vec[index+Q] <- stat.outputs$sd
    iv.str[index] <- stat.outputs$iv.str; iv.thol[index] <- stat.outputs$iv.thol;
    explained.iv[index] <- stat.outputs$explained.iv
    D.resid[[index]] <- stat.outputs$D.resid
    trace.T[index] = stat.outputs$trace.T
    RSS.vec[[index]] = stat.outputs$RSS.V
  }
  # Residual sum of squares of D.rep~Cov.rep
  D.RSS <- iv.str*SigmaSqD
  
  
  # violation space selection
  # all of the q below are from 0 to (Q-1), so use q+1 to index the columns
  # Comparison and robust estimators
  Coef.robust <- sd.robust <- rep(NA,4)
  names(Coef.robust) <- names(sd.robust) <- c("RF-comp","RF-Cor-comp","RF-robust","RF-Cor-robust")
  ivtest.vec <- (iv.str>=iv.thol)
  run.OLS <- weak.iv <- FALSE
  if (sum(ivtest.vec)==0) {
    warning("Weak IV, even if the IV is assumed to be valid; run OLS") # stop, output results of OLS
    run.OLS <- TRUE
    Qmax <- 1
  } else {
    Qmax <- sum(ivtest.vec)-1
    if (Qmax==0) {
      warning("Weak IV, if the IV is invalid. We still test the IV invalidity.") # zijian: we need to rewrite this sentence...
      Qmax <- 1
      weak.iv = TRUE
    }
  }
  # calculate Covariance at Qmax and compute bias correction
  eps.Qmax = eps.hat[[Qmax+1]]
  Coef.Qmax = rep(NA,Q)
  for (i in 1:Q) {
    Coef.Qmax[i] = Coef.vec[i] - sum(RSS.vec[[i]]*delta.hat*eps.Qmax)/D.RSS[i]
    Coef.vec[i+Q] = Coef.vec[i] - sum(RSS.vec[[i]]*delta.hat*eps.hat[[i]])/D.RSS[i]
  }
  sd.vec[-(1:Q)] = sd.vec[1:Q]
  
  ### Selection
  # define comparison matrix
  H <- beta.diff <- matrix(0,Qmax,Qmax)
  # compute H matrix
  for (q1 in 0:(Qmax-1)) {
    for (q2 in (q1+1):(Qmax)) {
      H[q1+1,q2] <- as.numeric(sum((weight%*%D.resid[[q1+1]])^2*eps.Qmax^2)/(D.RSS[q1+1]^2) + 
                                 sum((weight%*%D.resid[[q2+1]])^2*eps.Qmax^2)/(D.RSS[q2+1]^2) - 
                                 2*sum(eps.Qmax^2*(weight%*%D.resid[[q1+1]])*(weight%*%D.resid[[q2+1]]))/(D.RSS[q1+1]*D.RSS[q2+1])
      )
      
    }
  }
  # compute beta difference matrix, use Qmax
  for (q in 0:(Qmax-1)) {
    beta.diff[q+1,(q+1):(Qmax)] <- abs(Coef.Qmax[q+1]-Coef.Qmax[(q+2):(Qmax+1)]) # use bias-corrected estimator
  }
  # bootstrap for the quantile of the differences
  max.val <- rep(NA,300)
  eps.Qmax.cent = eps.Qmax - mean(eps.Qmax)
  for (i in 1:300) {
    diff.mat <- matrix(0,Qmax,Qmax)
    eps <- rep(NA,n.A1)
    for (j in 1:n.A1) {
      U.j = rnorm(1)
      eps[j] = eps.Qmax.cent[j]*U.j
    }
    eps.rep <- weight%*%eps
    for (q1 in 0:(Qmax-1)) {
      for (q2 in (q1+1):(Qmax)) {
        diff.mat[q1+1, q2] <- sum(D.resid[[q2+1]]*eps.rep)/(D.RSS[q2+1])-sum(D.resid[[q1+1]]*eps.rep)/(D.RSS[q1+1])
      }
    }
    diff.mat <- abs(diff.mat)/sqrt(H)
    max.val[i] <- max(diff.mat,na.rm = TRUE)
  }
  z.alpha <- 1.01*quantile(max.val,0.975)
  diff.thol <- z.alpha*sqrt(H)
  # comparison matrix
  C.alpha <- ifelse(beta.diff<=diff.thol,0,1)
  
  # layer selection or not
  if (layer==TRUE) {
    # a vector indicating the selection of each layer
    sel.vec <- apply(C.alpha,1,sum)
    if (all(sel.vec != 0)) {
      q.comp = Qmax
    } else {
      q.comp = min(which(sel.vec==0))-1
    }
  } else {
    ### What if Q = 3(q2) and Qmax = 1?
    sel.val <- C.alpha[1,Qmax]
    if (sel.val==1) {
      q.comp = Qmax
    } else {
      q.comp = 0
    }
  }
  
  ### invalidity of TSLS
  if (q.comp>=1) {
    invalidity <- 1
  } else {
    invalidity <- 0
  }
  q.robust <- min(q.comp+1, Qmax)
  Coef.robust[1] <- Coef.vec[q.comp+1]
  Coef.robust[2] <- Coef.vec[q.comp+Q+1]
  Coef.robust[3] <- Coef.vec[q.robust+1]
  Coef.robust[4] <- Coef.vec[q.robust+Q+1]
  sd.robust[1] <- sd.vec[q.comp+1]
  sd.robust[2] <- sd.vec[q.comp+Q+1]
  sd.robust[3] <- sd.vec[q.robust+1]
  sd.robust[4] <- sd.vec[q.robust+Q+1]
  CI.robust = rbind(Coef.robust + qnorm(alpha/2)*sd.robust,Coef.robust + qnorm(1-alpha/2)*sd.robust)
  rownames(CI.robust) = c("lower","upper")
  
  returnList = list(Coef.vec = Coef.vec,
                    sd.vec = sd.vec,
                    Coef.robust = Coef.robust,
                    sd.robust = sd.robust,
                    CI.robust = CI.robust,
                    iv.str = iv.str, iv.thol = iv.thol,
                    SigmaSqD = SigmaSqD,
                    SigmaSqY = SigmaSqY,
                    SigmaSqY.Qmax = mean(eps.Qmax^2),
                    trace.T = trace.T,
                    explained.iv = explained.iv,
                    Qmax = Qmax,
                    q.comp =q.comp, q.robust = q.robust,
                    invalidity = invalidity,
                    run.OLS = run.OLS,
                    weak.iv = weak.iv)
  returnList
}

