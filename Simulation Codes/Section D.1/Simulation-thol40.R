source("Source-RF-hetero2.R")

# parameter change -------------------------------------------------------

round = 1  # {1, 2, ..., 25}
# nonliearity degree
a = 0.35  # {0.15, 0.17, 0.19, 0.21, 0.23, 0.25, ..., 0.35}

myFile = file("stdin")
open(myFile)
round = as.numeric(readLines(myFile, n = 1))  
a = as.numeric(readLines(myFile, n = 1))  
print(c(
    round = round
    , a = a
))

# simulation --------------------------------------------------------------

sim_func = function(a) {
    # generate data -----------------------------------------------------------
    
    # dimension
    p = 20
    n = 3000
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
    reg.rf = lm(Y.rep ~ D.rep + Cov.rep)
    coef_tsci0 = coef(reg.rf)[2]
    eps.hat = resid(lm(Y[A1.ind,]-D[A1.ind,]*coef_tsci0 ~ Cov.aug[A1.ind,]))
    delta.hat = D[A1.ind,] - D.rep
    stat.outputs = TSCI.RF.stat(D.rep, Cov.rep, weight, n, eps.hat, delta.hat, str.thol = 10)
    ivstr = stat.outputs$iv.str
    RSS.vec = stat.outputs$RSS.V
    coef_tsci = coef_tsci0 - sum(RSS.vec*delta.hat*eps.hat) / (ivstr*mean(delta.hat^2))
    sd_tsci = stat.outputs$sd
    
    return(c(coef_tsci, sd_tsci, ivstr))
}

nsim = 20
total = NULL
for (i in 1:nsim) {
    print(i)
    result = sim_func(a)
    total = rbind(total, result)
}
filename = paste0("results-thol40/results",
                  "-a", a,
                  "-round", round, ".rds")
saveRDS(total, filename)



    

