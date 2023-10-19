library(RobustIV)
source("Source-RF-hetero2.R")
source("Source-CIIV.R")
source("Source-CIIV_helpers.R")

# parameter change --------------------------------------------------------

round = 1  # {1, 2, ..., 5}
# number of valid IVs
a = 0.5  # {0.1, 0.2, ..., 1}

# simulation --------------------------------------------------------------

sim_func <- function(a) {
    # true parameter
    beta0 = 1
    # sample size
    n = 3000
    # number of covariates
    p = 5
    # number of IVs
    L = 10
    # number of valid IVs
    l = 4
    # valid/invalid vector
    pi = a * c(rep(0,l), rep(0.5,L - l))
    # generate X and Z
    mu_Xstar = rep(0, p + L)
    cov_Xstar = stats::toeplitz(0.5^c(0:(p + L - 1)))
    Xstar = MASS::mvrnorm(n, mu_Xstar, cov_Xstar)
    Z = Xstar[, 1:L]
    X = Xstar[, (L + 1):(L + p)]
    # generate error
    Error_cov = matrix(c(1, 0.5, 0.5, 1), 2, 2)
    Error = MASS::mvrnorm(n, rep(0, 2), Error_cov)
    delta = Error[,1]
    epsilon = Error[,2]
    # generate D
    # D = rowSums(Z^2) + rowSums(tanh(Z)) + 0.5*rowSums(X) + delta  # all works, TSCI is more efficient and robust
    D = rowSums(abs(Z)) + rowSums(tanh(Z)) + 0.5*rowSums(X) + delta  # all works, they are more efficient, TSCI is more robust
    # generate Y
    Y = D*beta0 + Z %*% pi + 0.5*rowSums(X) + epsilon
    
    # TSHT --------------------------------------------------------------------
    
    res_tsht = TSHT(Y, D, Z, X)
    
    # CIIV --------------------------------------------------------------------
    
    CIIV_temp = tryCatch_E(CIIV(Y,D,Z,X,robust = FALSE,firststage = T), 
                           ret.obj = list(NULL))
    if (is.null(CIIV_temp$error)) {
        res_CIIV = CIIV_temp$value
    }
    
    # TSCI --------------------------------------------------------------------
    
    Y <- as.matrix(Y)
    D <- as.matrix(D)
    Z <- as.matrix(Z)
    X <- as.matrix(X)
    n <- nrow(X)
    P <- ncol(X) + ncol(Z)
    
    V1 = Z
    V2 = Z^2
    vio.space = list(V2, V1)
    
    res_tsci = TSCI.RF(Y,D,Z,X
                       ,vio.space = vio.space
                       ,mtry = seq(round(P/3),round(2*P/3),by = 1)
                       ,min.node.size = c(5,10,20)
                       ,num.trees = 500)
    
    # return ------------------------------------------------------------------
    
    returnResults = list(res_tsht = res_tsht,
                         res_CIIV = res_CIIV,
                         res_tsci = res_tsci)
    returnResults
}

nsim = 20
total = NULL
for (i in 1:nsim) {
    print(i)
    result = sim_func(a)
    total = rbind(total, result)
}
filename = paste0("results-noRule/results",
                  "-a", a,
                  "-round", round, ".rds")
saveRDS(total, filename)

