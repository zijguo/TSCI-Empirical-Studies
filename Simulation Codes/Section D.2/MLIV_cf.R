# get the coefficient of MLIV in the D~fHat+X
source('Source-RF-hetero2.R')
library(ranger)
library(Matrix)
library(MASS)
library(ivreg)

sim_func <- function(n, a) {
    # generate data -----------------------------------------------------------
    
    # dimension
    p = 20
    # X and Z
    mu_Xstar = rep(0, p + 1)
    cov_Xstar = stats::toeplitz(0.5^c(0:p))
    Xstar = MASS::mvrnorm(n, mu_Xstar, cov_Xstar)
    X = pnorm(Xstar)[,-1]
    # Z = (pnorm(Xstar)[,1] - 0.5) * 4  # ivstr too large
    Z = pnorm(Xstar)[,1]
    # error
    delta = epsilon = matrix(NA, nrow = n, ncol = 1)
    tau1 = rep(0, n)
    for (j in 1:n) {
        tau1[j] = rnorm(1, 0, sqrt(Z[j]^2 + 0.25))
    }
    tau2 = rnorm(n)
    for (j in 1:n) {
        delta[j,1] = rnorm(1, 0, sqrt(Z[j]^2 + 0.25))
        epsilon[j,1] = 0.6*delta[j,1] + sqrt((1 - 0.6^2) / (0.86^4 + 1.38072^2)) * (tau1[j]*1.38072 + tau2[j]*0.86^2) 
    }
    # D
    f_func = function(Z, X, a) {
        0.5*Z + a*(sin(2*pi*Z) + 1.5*cos(2*pi*Z)) - 0.3*rowSums(X[,1:10])
    }
    D = f_func(Z, X, a) + delta
    # Y
    phi_func = function(X) 0.2*rowSums(X[,1:10])
    vio = 0
    if (vio == 0) {
        h_func = function(Z, X) 0
    } else if (vio == 1) {
        h_func = function(Z, X) Z
    } else if (vio == 2) {
        h_func = function(Z, X) Z + Z^2 - 1
    }
    Y = D*0.5 + h_func(Z, X) + phi_func(X) + epsilon
    
    # TSCI --------------------------------------------------------------------
    
    Y = as.matrix(Y)
    D = as.matrix(D)
    Z = as.matrix(Z)
    X = as.matrix(X)
    n = nrow(X)
    P = ncol(X) + ncol(Z)
    Cov.aug = X
    n.A1 = round(2/3*n)
    A1.ind = 1:n.A1
    forest = TSCI.RF.fit(D,Z,X,
                         num.trees = 200,
                         mtry = seq(round(P/3),round(2*P/3),by = 1),
                         max.depth = 0,
                         min.node.size = c(5,10,20),
                         split.prop = 2/3)
    weight = TSCI.RF.weight(forest$nodes.A1)
    Y.rep = as.matrix(weight %*% Y[A1.ind,])
    D.rep = as.matrix(weight %*% D[A1.ind,])
    Cov.rep = as.matrix(weight %*% Cov.aug[A1.ind,])
    
    # estimator decomposition
    tsci_deno = sum(resid(lm(D.rep ~ Cov.rep))^2)
    tsci_nume = sum(resid(lm(D.rep ~ Cov.rep)) * resid(lm(Y.rep ~ Cov.rep)))
    coef_tsci0 = tsci_nume / tsci_deno

    eps.hat = resid(lm(Y[A1.ind,]-D[A1.ind,]*coef_tsci0 ~ Cov.aug[A1.ind,]))
    delta.hat = D[A1.ind,] - D.rep
    stat.outputs = TSCI.RF.stat(D.rep, Cov.rep, weight, n, eps.hat, delta.hat, str.thol = 10)
    ivstr = stat.outputs$iv.str
    RSS.vec = stat.outputs$RSS.V
    coef_tsci = coef_tsci0 - sum(RSS.vec*delta.hat*eps.hat) / (ivstr*mean(delta.hat^2))
    
    # MLIV --------------------------------------------------------------------
    
    # res_mliv = ivmodel(Y[A1.ind,], D[A1.ind,], D.rep, X[A1.ind,], heteroSE = T)
    # coef_mliv = coef(res_mliv)["TSLS", "Estimate"]
    Y.A1 = as.matrix(Y[A1.ind,])
    D.A1 = as.matrix(D[A1.ind,])
    Z.A1 = as.matrix(Z[A1.ind,])
    X.A1 = as.matrix(X[A1.ind,])
    Cov.aug.A1 = as.matrix(Cov.aug[A1.ind,])
    stage1 = lm(D.A1 ~ D.rep + X.A1)
    stage2 = lm(Y.A1 ~ predict(stage1) + Cov.aug.A1)
    
    mliv_deno = sum(resid(lm(predict(stage1) ~ Cov.aug.A1))^2)
    c_f = coef(stage1)[2]
    mliv_deno_nocf = sum(resid(lm(D.rep ~ Cov.aug.A1))^2)
    
    mliv_nume = sum(resid(lm(predict(stage1) ~ Cov.aug.A1)) * resid(lm(Y.A1 ~ Cov.aug.A1)))
    mliv_nume_nocf = sum(resid(lm(D.rep ~ Cov.aug.A1)) * resid(lm(Y.A1 ~ Cov.aug.A1)))
    
    # return ------------------------------------------------------------------
    
    returnValue = list(
        tsci_deno = tsci_deno,
        tsci_nume = tsci_nume,
        coef_tsci = coef_tsci,
        ivstr = ivstr,
        mliv_deno = mliv_deno,
        mliv_nume = mliv_nume,
        c_f = c_f,
        mliv_deno_nocf = mliv_deno_nocf,
        mliv_nume_nocf = mliv_nume_nocf
    )
    returnValue

}


n = 1000
a = 0.2
total = NULL
for (i in 1:500) {
    print(i)
    cf = sim_func(n, a)
    total = rbind(total, cf)
}
saveRDS(total, "results-MLIV_cf.rds")

