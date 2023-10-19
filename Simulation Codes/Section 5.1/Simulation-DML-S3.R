library(ivmodel)
library(fda)
library(DoubleML)
library(mlr3)
library(mlr3learners)
library(mlr3tuning)
library(data.table)
source("Source-DML_helpers.R")
source("Source-RF-hetero2.R")

# parameter change --------------------------------------------------------

round = 1  # {1, 2, ..., 25}
# concen
a = 0.4  # {0.2, 0.4, ..., 1.6}

# simulation --------------------------------------------------------------

sim_func = function(a) {
    ## Generate some data
    # dimension
    p = 5
    n = 3000
    # X and Z
    # mu_Xstar = rep(0, p+1)
    # cov_Xstar = stats::toeplitz(0.5^c(0:p))
    # Xstar = MASS::mvrnorm(n, mu_Xstar, cov_Xstar)
    # X = pnorm(Xstar)[,-1]
    X = matrix(runif(n*p, -1, 1), ncol = p)
    Z = 3*tanh(2*X[,1]) + rnorm(n, 0, 1)
    # hidden confounders
    h <- 2*sin(X[,1]) + rnorm(n, 0, 1)
    # D
    f_func = function(Z, X, a) {
        # Z + a*(sin(2*pi*Z) + 1.5*cos(2*pi*Z) + 0.1*Z*rowSums(X[,1:10])) - 0.5*rowSums(X)  # model1: seems good
        # Z/2 + a*(sin(2*pi*Z)+1.5*cos(2*pi*Z)+Z^3/8) - 0.2*rowSums(X)                      # model2: not good, ivstr<concen
        a*Z + sin(2*pi*Z) + 1.5*cos(2*pi*Z) - 0.2*rowSums(X)                              # model3: a up to 2, good, ivstr always>1000, concen from 20 to 1000
        # a*(Z/4+Z*rowSums(X[,1:5])/2) + sin(2*pi*Z) + cos(2*pi*Z) - 0.2*rowSums(X)         # model4: a up to 0.7
        
    }
    D = f_func(Z, X, a) - 0.5 * h + rnorm(n, 0, 1)
    # Y
    phi_func = function(X) X[,1]*X[,2] + X[,2]*X[,3] + X[,3]*X[,4] + X[,4]*X[,5] + X[,5]*X[,1]
    Y = D*0.5 + phi_func(X) - cos(pi * 0.25 * h) + rnorm(n, 0, 1)
    
    # DML ---------------------------------------------------------------------
    
    ml_l = lrn("regr.ranger", num.trees = 500, max.depth = 0)
    ml_m = ml_l$clone()
    ml_r = ml_l$clone()
    dt = data.table(cbind(Y, D, Z, X))
    colnames(dt) = c("y","d","z",paste0("x",1:ncol(X)))
    obj_dml_data = DoubleMLData$new(dt, y_col = "y", d_cols = "d", 
                                    x_cols = paste0("x",1:ncol(X)), z_cols = "z")
    dml_pliv_obj = DoubleMLPLIV$new(obj_dml_data, ml_l, ml_m, ml_r)
    P = ncol(as.matrix(X)) + 1
    param_grid = list(
        "ml_l" = paradox::ParamSet$new(list(
            paradox::ParamInt$new("mtry", lower = round(round(P/3)), upper = round(2*P/3)),
            paradox::ParamInt$new("min.node.size", lower = 5, upper = 20, special_vals = list(5,10,20)))),
        "ml_m" = paradox::ParamSet$new(list(
            paradox::ParamInt$new("mtry", lower = round(round(P/3)), upper = round(2*P/3)),
            paradox::ParamInt$new("min.node.size", lower = 5, upper = 20, special_vals = list(5,10,20)))),
        "ml_r" = paradox::ParamSet$new(list(
            paradox::ParamInt$new("mtry", lower = round(round(P/3)), upper = round(2*P/3)),
            paradox::ParamInt$new("min.node.size", lower = 5, upper = 20, special_vals = list(5,10,20)))))
    # paradox::ParamInt$new("min.node.size", lower = 1, upper = 2))))
    tune_settings = list(
        terminator = mlr3tuning::trm("evals", n_evals = 5),
        algorithm = mlr3tuning::tnr("grid_search", resolution = 5))
    dml_pliv_obj$tune(param_set = param_grid, tune_settings = tune_settings)
    dml_pliv_obj$fit()
    beta_DMLorg = dml_pliv_obj$coef
    sd_DMLorg = dml_pliv_obj$se
    
    # TSLS --------------------------------------------------------------------
    
    # same as the above, we use 2 decomposition ways
    model_tsls = ivmodel(Y, D, Z, X, heteroSE = T)
    beta_tsls = coef(model_tsls)["TSLS", "Estimate"]
    
    ### variance estimator
    eps_hat = resid(lm(Y-D*beta_tsls~X))
    sd_tsls = coef(model_tsls)["TSLS", "Std. Error"]
    
    ### concentration parameter
    stage1 = lm(D~Z+X)
    delta_hat = resid(stage1)
    nume = sum(resid(lm(predict(stage1) ~ X))^2)
    concen_tsls = nume / mean(delta_hat^2)
    
    # TSCI --------------------------------------------------------------------
    
    Y <- as.matrix(Y)
    D <- as.matrix(D)
    Z <- as.matrix(Z)
    X <- as.matrix(X)
    n <- nrow(X)
    P <- ncol(X) + ncol(Z)
    getDesign <- function(X, knots=5) {
        p <- NCOL(X)
        m = matrix(NA,nrow(X),0)
        for (j in 1:p) {
            X.current <- X[,j]
            UX.current <- unique(X.current)
            knots.use <- quantile(UX.current, seq(0, 1, length = knots))
            stopifnot(all.equal(range(knots.use), range(X.current)))
            basis <- create.bspline.basis(rangeval = range(knots.use),
                                          breaks = knots.use, norder = 4)
            m.current <- eval.basis(X.current, basis)
            m = cbind(m,m.current)
        }
        return(m)
    }
    Cov.aug <- getDesign(X)
    forest <- TSCI.RF.fit(D, Z, X,
                          num.trees = 500,
                          mtry = seq(round(P/3), round(2*P/3), by=1),
                          min.node.size = c(5,10,20),
                          max.depth = 0,
                          split.prop = 2/3)
    weight <- TSCI.RF.weight(forest$nodes.A1)
    A1.ind = forest$A1.ind
    Y.rep <- as.matrix(weight %*% Y[A1.ind])
    D.rep <- as.matrix(weight %*% D[A1.ind])
    Cov.rep <- as.matrix(weight %*% Cov.aug[A1.ind,])
    reg.rf <- lm(Y.rep ~ D.rep + Cov.rep)
    coef_tsci0 <- coef(reg.rf)[2]
    eps.hat = resid(lm(Y[A1.ind]-D[A1.ind]*coef_tsci0 ~ Cov.aug[A1.ind,]))
    delta.hat = D[A1.ind] - D.rep
    stat.outputs <- TSCI.RF.stat(D.rep,Cov.rep,weight,n,eps.hat,delta.hat,str.thol = 10)
    RSS.vec <- stat.outputs$RSS.V
    ivstr = stat.outputs$iv.str
    beta_tsci = coef_tsci0 - sum(RSS.vec*delta.hat*eps.hat) / (ivstr*mean(delta.hat^2))
    
    ### variance estimator
    sd_tsci = stat.outputs$sd
    
    returnRecord = list(
        # DMLorg
        beta_DMLorg = beta_DMLorg,
        sd_DMLorg = sd_DMLorg,
        # TSLS
        beta_tsls = beta_tsls,
        sd_tsls = sd_tsls,
        concen_tsls = concen_tsls,
        # TSCI
        beta_tsci = beta_tsci,
        beta_tsci0 = coef_tsci0,
        sd_tsci = sd_tsci,
        ivstr = ivstr
    )
    returnRecord
}

nsim = 20
total = NULL
for (i in 1:nsim) {
    print(i)
    result = sim_func(a)
    total = rbind(total, result)
}
filename = paste0("results-misspecified-TSCIworks/results",
                  "-a", a,
                  "-round", round, ".rds")
saveRDS(total, filename)


