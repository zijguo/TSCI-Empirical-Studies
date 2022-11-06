# use f hat as IV in TSLS and CF and EE (valid IV setting)
source('Source-RF-hetero2.R')
library(ranger)
library(Matrix)
library(MASS)
library(ivreg)

nsim <- 100

sim_func_fhat <- function(n, a) {
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
  weight <- TSCI.RF.weight(forest$nodes.A1)
  
  
  # TSLS --------------------------------------------------------------------
  
  Y_A1 <- as.matrix(Y[A1_ind])
  D_A1 <- as.matrix(D[A1_ind])
  X_A1 <- as.matrix(X[A1_ind,])
  f_A1 <- as.matrix(f[A1_ind])
  Z_A1 <- as.matrix(Z[A1_ind])
  f_hat <- as.matrix(weight %*% D_A1)
  
  reg <- ivreg(Y_A1 ~ D_A1 + X_A1 | f_hat + X_A1)
  beta_TSLS <- coef(reg)[2]
  sd_TSLS <- sqrt(diag(vcov(reg))[2])
  CI_TSLS <- confint(reg)[2, ]
  cover_TSLS <- as.numeric(CI_TSLS[1] <= beta & beta <= CI_TSLS[2])
  
  concen_param <- function(Z, X, ITT_D, SigmaSqD){
    qrX = qr(cbind(X,1))  # the input X has no constant column, and we do Z~X with intercept
    resid_Z = qr.resid(qrX, Z)
    numerator = sum((resid_Z %*% ITT_D)^2)
    return(numerator / SigmaSqD)
  }
  stage1 <- lm(D_A1 ~ f_hat + X_A1)
  SigmaSqD <- mean(resid(stage1)^2)
  concen <- concen_param(f_hat, X_A1, as.matrix(coef(stage1)[2]), SigmaSqD)
  
  
  
  # Control function --------------------------------------------------------
  
  stage1 <- lm(D_A1 ~ f_hat + X_A1)
  stage2 <- lm(Y_A1 ~ D_A1 + X_A1 + resid(stage1))
  beta_CF <- coef(stage2)[2]
  eps_hat <- resid(lm(Y_A1-D_A1*beta_CF ~ X_A1))
  D_pro <- resid(lm(D_A1 ~ X_A1 + resid(stage1)))
  sd_CF <- sqrt(solve(crossprod(D_pro)) * mean(eps_hat^2))
  CI_CF <- c(beta_CF+qnorm(0.05/2)*sd_CF, beta_CF+qnorm(1-0.05/2)*sd_CF)
  cover_CF <- as.numeric(CI_CF[1] <= beta & beta <= CI_CF[2])
  
  
  # Moment equation (betaEE) ---------------------------------------------------------
  
  Cov_A1 <- as.matrix(X[A1_ind, ])
  f_hat_resid <- resid(lm(f_hat ~ Cov_A1))
  Y_resid <- resid(lm(Y_A1 ~ Cov_A1))
  D_resid <- resid(lm(D_A1 ~ Cov_A1))
  beta_EE <- sum(Y_resid * f_hat_resid) / sum(D_resid * f_hat_resid)
  eps_hat <- resid(lm(Y_A1 - D_A1*beta_EE ~ Cov_A1))
  sd_EE <- sqrt(sum(f_hat_resid^2) * mean(eps_hat^2) / sum(D_resid * f_hat_resid)^2)
  CI_EE <- c(beta_EE+qnorm(0.05/2)*sd_EE, beta_EE+qnorm(1-0.05/2)*sd_EE)
  cover_EE <- as.numeric(CI_EE[1] <= beta & beta <= CI_EE[2])
  
  return(c(
    beta_TSLS = as.numeric(beta_TSLS)
    , bias_TSLS = as.numeric(beta_TSLS - beta)
    , sd_TSLS = as.numeric(sd_TSLS)
    , cover_TSLS = cover_TSLS
    , length_TSLS = as.numeric(CI_TSLS[2] - CI_TSLS[1])
    , concen = concen
    , beta_CF = as.numeric(beta_CF)
    , bias_CF = as.numeric(beta_CF - beta)
    , sd_CF = as.numeric(sd_CF)
    , cover_CF = cover_CF
    , length_CF = as.numeric(CI_CF[2] - CI_CF[1])
    , beta_EE = beta_EE
    , bias_EE = beta_EE - beta
    , sd_EE = sd_EE
    , cover_EE = cover_EE
    , length_EE = CI_EE[2] - CI_EE[1]
  ))
}

total <- NULL
for (a in c(0, 0.1, 0.125, 0.15, 0.2, 0.5, 0.75, 1)) {
  for (n in c(1000, 3000, 5000)) {
    records <- NULL
    for (i in seq(nsim)) {
      print(c(a, n, i))
      result <- sim_func_fhat(n, a)
      records <- rbind(records, result)
    }
    total <- rbind(total, c(a=a, n=n, colMeans(records), sd_emp_TSLS=sd(records[,1]), sd_emp_CF=sd(records[,7]), sd_emp_EE=sd(records[,12])))
  }
}

