source("Source-RF-hetero2.R")

# parameter change -------------------------------------------------------

round = 1  # {1, 2, ..., 25}
# sample size
n = 1000  # {1000, 3000, 5000}
# nonliearity degree
a = 0.25  # {0.25, 0.5, 0.75}

# simulation --------------------------------------------------------------

sim_func <- function(n, a) {
    # generate data -----------------------------------------------------------
    
    # dimension
    p = 5
    # X and Z
    mu_Xstar = rep(0, p+1)
    cov_Xstar = stats::toeplitz(0.5^c(0:p))
    Xstar = MASS::mvrnorm(n, mu_Xstar, cov_Xstar)
    X = Xstar[,-1]
    Z = pnorm(Xstar)[,1]
    Z = ifelse(Z > 0.6, 1, 0)
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
        # -25/12 + Z + Z^3/3 + a*Z*rowSums(X[,1:5]) - 0.3*rowSums(X)
        Z + a*(Z*rowSums(X[,c(1,2,3,4)]) + Z*rowSums(X[,c(1,2,3,4)]*X[,c(2,3,4,5)])) - 0.3*rowSums(X[,1:5])
    }
    D = f_func(Z, X, a) + delta
    # Y
    h_func = function(Z, X) {
        Z + 0.5*Z*rowSums(X[,1:3])
    }
    Y = D + h_func(Z, X) + 0.2*rowSums(X) + epsilon
    
    # TSCI --------------------------------------------------------------------
    
    Q = 2
    Z = as.matrix(Z)
    interaction = matrix(NA,n,5)
    for (j in 1:5) {
        interaction[,j] = Z*X[,j]
    }
    vio.space = list(cbind(interaction,Z))
    
    res_tsci = TSCI.RF(Y,D,Z,X,
                       vio.space = vio.space,
                       mtry=c(1:p),
                       min.node.size=c(1,3,5,8,10,12,14,16,18,20))
    
    # return ------------------------------------------------------------------
    
    res_tsci
}

nsim = 20
total = NULL
for (i in 1:nsim) {
    print(i)
    result = sim_func(n, a)
    total = rbind(total, result)
}
filename = paste0("results-model3/results-Simulation-TSCI-invalidIV-model3"
                  , "-a", a
                  , "-n", n
                  , "-round", round, ".rds")
saveRDS(total, filename)

