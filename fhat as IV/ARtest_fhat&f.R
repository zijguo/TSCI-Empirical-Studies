# Anderson-Rubin test: use fhat and f as IV
source('Source-RF-hetero2.R')
library(ranger)
library(Matrix)
library(MASS)
library(ivmodel)


nsim <- 500
ar_sim <- function(n, a) {
  # generate data -----------------------------------------------------------
  
  # dimension
  p <- 10
  # sample size
  n <- n
  # interaction value
  inter_val <- 0
  # the IV strength
  a <- a
  # violation strength
  tau <- 0
  f_1 <- function(z, a) {
    0.5*z + a * ( z^2 + 0.125 * z^4) - 25 / 12
  }
  rho <- 0.5
  Cov <- stats::toeplitz(rho^c(0 : p))
  mu <- rep(0, p + 1)
  # true effect
  beta <- 1
  alpha <- as.matrix(rep(-0.3, p))
  gamma <- as.matrix(rep(0.2, p))
  inter <- as.matrix(c(rep(inter_val, 5),rep(0, p - 5)))
  
  
  # generate the data
  mu_error <- rep(0,2)
  Cov_error <- matrix(c(1, 0.5, 0.5,1), 2, 2)
  Error <- MASS::mvrnorm(n, mu_error, Cov_error)
  W_original <- MASS::mvrnorm(n, mu, Cov)
  W <- pnorm(W_original)
  # instrument variable
  Z <- W[, 1]
  # baseline covariates
  X <- W[, -1]
  # generate the treatment variable D
  D <- f_1(Z, a) + X %*% alpha + Z * X %*% inter + Error[, 1]
  f <- D - Error[,1]
  # generate the outcome variable Y
  Y <- D * beta + tau * Z + X %*% gamma + Error[, 2]
  
  
  
  # fitting f hat -----------------------------------------------------------
  
  forest <- TSCI.RF.fit(D,Z,X
                        ,num.trees=200
                        ,mtry=seq(round((ncol(X)+1)/3),round(2*(ncol(X)+1)/3),by=1)
                        ,max.depth=0
                        ,min.node.size=c(5,10,20)
                        ,split.prop=2/3)
  A1_ind <- forest$A1.ind
  # Data_A1 <- as.data.frame(cbind(D, Z, X))[A1_ind, ]
  # names(Data_A1) <- c("D", paste("W", 1:(ncol(X)+1), sep = ""))
  # f_hat <- predict(forest$forest.A2, data = Data_A1, type = "response")$predictions
  TSCI.RF.weight.wtself <- function(nodes) {
    n.A1 <- NROW(nodes); num.trees <- NCOL(nodes)
    out.weight <- matrix(0,n.A1,n.A1)
    for (j in 1:num.trees) {
      weight.mat <- matrix(0,n.A1,n.A1) # weight matrix for single tree
      unique.nodes <- unique(nodes[,j])
      for (i in 1:length(unique.nodes)) {
        ind <- nodes[,j]==unique.nodes[i] # indices of samples in the node
        num.samples <- sum(ind) # number of samples in the node
        w <- 1/(num.samples)  # weight, to remove self-prediction
        weight.vec <- ifelse(ind,yes=w,no=0)
        weight.mat[ind,] <- matrix(rep(weight.vec,num.samples),num.samples,byrow=T)/num.trees
      }
      # diag(weight.mat) <- 0 # remove self prediction
      out.weight <- out.weight + weight.mat
    }
    out.weight <- Matrix(out.weight, sparse = T) # sparse matrix to save memory
    return(out.weight)
  }
  
  weight_no <- TSCI.RF.weight(forest$nodes.A1)
  weight_wt <- TSCI.RF.weight.wtself(forest$nodes.A1)
  
  
  
  # TSLS and AR test -----------------------------------------------------------------
  
  Y_A1 <- as.matrix(Y[A1_ind])
  D_A1 <- as.matrix(D[A1_ind])
  X_A1 <- as.matrix(X[A1_ind,])
  f_A1 <- as.matrix(f[A1_ind])
  Z_A1 <- as.matrix(Z[A1_ind])
  f_hat_no <- as.matrix(weight_no %*% D_A1)
  f_hat_wt <- as.matrix(weight_wt %*% D_A1)
  
  # fhat_no
  reg <- ivmodel(Y_A1, D_A1, f_hat_no, X_A1, k=1)
  beta_fhat_no <- coef(reg)['TSLS', 'Estimate']
  sd_fhat_no <- coef(reg)['TSLS', 'Std. Error']
  CI_fhat_no <- reg$AR$ci
  length_fhat_no <- sum(CI_fhat_no[, 2] - CI_fhat_no[, 1])
  cover_fhat_no <- any(CI_fhat_no[, 1] <= beta & beta <= CI_fhat_no[, 2])
  
  # fhat_wt
  reg <- ivmodel(Y_A1, D_A1, f_hat_wt, X_A1, k=1)
  beta_fhat_wt <- coef(reg)['TSLS', 'Estimate']
  sd_fhat_wt <- coef(reg)['TSLS', 'Std. Error']
  CI_fhat_wt <- reg$AR$ci
  length_fhat_wt <- sum(CI_fhat_wt[, 2] - CI_fhat_wt[, 1])
  cover_fhat_wt <- any(CI_fhat_wt[, 1] <= beta & beta <= CI_fhat_wt[, 2])
  
  # f
  reg <- ivmodel(Y_A1, D_A1, f_A1, X_A1, k=1)
  beta_f <- coef(reg)['TSLS', 'Estimate']
  sd_f <- coef(reg)['TSLS', 'Std. Error']
  CI_f <- reg$AR$ci
  length_f <- sum(CI_f[, 2] - CI_f[, 1])
  cover_f <- any(CI_f[, 1] <= beta & beta <= CI_f[, 2])
  
  # concentration parameter
  concen_param <- function(Z, X, ITT_D, SigmaSqD){
    qrX = qr(cbind(X,1))  # the input X has no constant column, and we do Z~X with intercept
    resid_Z = qr.resid(qrX, Z)
    numerator = sum((resid_Z %*% ITT_D)^2)
    return(numerator / SigmaSqD)
  }
  stage1 <- lm(D_A1 ~ f_hat_no + X_A1)
  SigmaSqD <- mean(resid(stage1)^2)
  concen_no <- concen_param(f_hat_no, X_A1, as.matrix(coef(stage1)[2]), SigmaSqD)
  
  stage1 <- lm(D_A1 ~ f_hat_wt + X_A1)
  SigmaSqD <- mean(resid(stage1)^2)
  concen_wt <- concen_param(f_hat_wt, X_A1, as.matrix(coef(stage1)[2]), SigmaSqD)
  
  
  return(c(
    bias_fhat_no = beta_fhat_no - beta
    , bias_fhat_wt = beta_fhat_wt - beta
    , bias_f = beta_f - beta
    , sd_fhat_no = sd_fhat_no
    , sd_fhat_wt = sd_fhat_wt
    , sd_f = sd_f
    , cover_fhat_no = cover_fhat_no
    , cover_fhat_wt = cover_fhat_wt
    , cover_f = cover_f
    , length_fhat_no = length_fhat_no
    , length_fhat_wt = length_fhat_wt
    , length_f = length_f
    , concen_no = concen_no
    , concen_wt = concen_wt
  ))
}

a <- 0
n <- 1000

total <- NULL
reshuffle_num <- 0
omit_num <- 0
for (i in seq(nsim)) {
  print(c(a, n, i))
  records <- try(ar_sim(n, a))
  if (inherits(records, "try-error")==T) {
    print('reshuffle')
    reshuffle_num <- reshuffle_num + 1
    records <- try(ar_sim(n, a))
    if (inherits(records, "try-error")==T) {
      print('omit')
      omit_num <- omit_num + 1
      next
    }
  }
  total <- rbind(total, records)
}


