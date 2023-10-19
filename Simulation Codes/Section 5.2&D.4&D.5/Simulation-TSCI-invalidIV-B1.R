library(ivmodel)
source("Source-otherRF-hetero.R")
source("Source-RF-hetero2.R")
# source("Source-otherRF-homo.R")

# parameter change -------------------------------------------------------

round = 1  # {1, 2, ..., 25}
# sample size
n = 1000  # {1000, 3000, 5000}
# nonliearity degree
a = 1  # {0, 0.5, 1}
# outcome model violation
vio = 1  # {1, 2}

# simulation --------------------------------------------------------------

sim_func <- function(n, a, vio) {
    # generate data -----------------------------------------------------------
    
    # dimension
    p = 20
    # X and Z
    mu_Xstar = rep(0, p+1)
    cov_Xstar = stats::toeplitz(0.5^c(0:p))
    Xstar = MASS::mvrnorm(n, mu_Xstar, cov_Xstar)
    X = pnorm(Xstar)[,-1]
    Z = (pnorm(Xstar)[,1] - 0.5) * 4
    # error
    delta = epsilon = matrix(NA, nrow=n, ncol=1)
    tau1 = rep(0, n)
    for (j in 1:n) {
        tau1[j] = rnorm(1, 0, sqrt(Z[j]^2+0.25))
    }
    tau2 = rnorm(n)
    for (j in 1:n) {
        delta[j,1] = rnorm(1, 0, sqrt(Z[j]^2+0.25))
        epsilon[j,1] = 0.6*delta[j,1] + sqrt((1-0.6^2)/(0.86^4+1.38072^2))*(tau1[j]*1.38072+tau2[j]*0.86^2) 
    }
    # D
    f_func = function(Z, X, a) {
        -25/12 + Z + Z^3/3 + a*Z*rowSums(X[,1:5]) - 0.3*rowSums(X)
    }
    D = f_func(Z, X, a) + delta
    # Y
    if (vio == 0) {
        h_func = function(Z, X) 0
    } else if (vio == 1) {
        h_func = function(Z, X) Z
    } else if (vio == 2) {
        h_func = function(Z, X) Z + Z^2 - 1
    }
    Y = D + h_func(Z, X) + 0.2*rowSums(X) + epsilon
    
    # TSCI --------------------------------------------------------------------
    
    res_tsci = TSCI.RF(Y, D, Z, X)
    
    # TSLS --------------------------------------------------------------------
    
    res_tsls = ivmodel(Y, D, Z, X, heteroSE=T)
    # concentration parameter
    stage1 = lm(D ~ Z + X)
    delta_hat = resid(stage1)
    nume = sum(resid(lm(predict(stage1) ~ X))^2)
    concen = nume / mean(delta_hat^2)

    # TSCI-plug ---------------------------------------------------------------
    
    Cov.aug = X
    for (k in 1:3) {
        Cov.aug = cbind(Z^k, Cov.aug)
    }
    forest.plug = TSCI.RF.fit(D, Z, X
                              ,num.trees=200
                              ,mtry=seq(round((p+1)/3), round(2*(p+1)/3), by=1)
                              ,max.depth=0
                              ,min.node.size=c(5,10,20)
                              ,split.prop=2/3)
    weight.plug = TSCI.RF.weight(forest.plug$nodes.A1)
    res_plug_hetero = naiveRF.stat.hetero(Y, D, Cov.aug[,-(1:(3-vio))], forest.plug$A1.ind, weight.plug)
    # res_plug = naiveRF.stat(Y, D, Cov.aug[,-(1:(3-vio))], forest.plug$A1.ind, weight.plug)
    rm(forest.plug, weight.plug)
    gc()
    
    # TSCI-full ---------------------------------------------------------------
    
    forest.full = TSRF.full(cbind(Z, X), D
                            ,mtry=seq(round((p+1)/3), round(2*(p+1)/3), by=1)
                            ,max.depth=0
                            ,min.node.size=c(5,10,20))
    weight.full = TSCI.RF.weight(forest.full$nodes)
    res_full_hetero = TSRF.stat.full.hetero(Y, D, Cov.aug[,-(1:(3-vio))], weight.full)
    # res_full = TSRF.stat.full(Y, D, Cov.aug[,-(1:(3-vio))], weight.full)
    rm(forest.full, weight.full)
    gc()

    # return ------------------------------------------------------------------
    
    returnValue = list(res_tsci = res_tsci
                       # , res_full = res_full
                       , res_full_hetero = res_full_hetero
                       # , res_plug = res_plug
                       , res_plug_hetero = res_plug_hetero
                       , res_tsls = c(coef(res_tsls)["TSLS",2], confint(res_tsls)["TSLS",], concen)
    )
    returnValue
}

nsim = 20
total = NULL
for (i in 1:nsim) {
    print(i)
    result = sim_func(n, a, vio)
    total = rbind(total, result)
}
filename = paste0("results-model1/results-Simulation-TSCI-invalidIV-model1"
                  , "-a", a
                  , "-vio", vio
                  , "-n", n
                  , "-round", round, ".rds")
saveRDS(total, filename)

