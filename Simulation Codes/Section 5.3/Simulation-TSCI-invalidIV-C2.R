library(ivmodel)
library(DoubleML)
library(mlr3)
library(mlr3learners)
library(mlr3tuning)
library(data.table)
source("Source-Basis-hetero.R")
source("Source-RF-hetero2.R")

# parameter change -------------------------------------------------------

round = 1  # {1, 2, ..., 25}
# nonliearity degree
a = 1.2  # {0, 0.1, 0.2, ..., 1.2}

# simulation --------------------------------------------------------------

sim_func <- function(a) {
    # generate data -----------------------------------------------------------
    
    # dimension
    p = 5
    n = 3000
    # X and Z
    X = matrix(pi*runif(n*p, -1, 1), ncol = p)
    Z = 3*tanh(2*X[,1]) + rnorm(n, 0, 1)
    # hidden confounders
    h <- 2*sin(X[,1]) + rnorm(n, 0, 1)
    # D
    f_func = function(Z, X, a) {
        -25/12 + Z + a*Z*rowSums(X[,1:5])/2 - 0.3*rowSums(X)
    }
    D = f_func(Z, X, a) - 0.5 * h + rnorm(n, 0, 1)
    # Y
    phi_func = function(X) rowSums(X^2)
    h_func = function(Z, X) Z
    Y = D*0.5 + h_func(Z, X) + phi_func(X) - cos(pi * 0.25 * h) + rnorm(n, 0, 1)
    
    # TSCI --------------------------------------------------------------------
    
    Y <- as.matrix(Y)
    D <- as.matrix(D)
    Z <- as.matrix(Z)
    X <- as.matrix(X)
    P <- ncol(X) + ncol(Z)
    # define the vio.space as polynomials if not specified
    Q = 4
    vio.space <- matrix(NA,n,0)
    for (q in 1:(Q - 1)) {
        vio.space <- cbind(Z^q,vio.space)
    }
    # the indices to remove to identify violation space
    rm.ind = rep(list(NA),Q - 1)
    for (i in 1:(Q - 1)) {
        rm.ind[[i]] = 1:(Q - i)
    }
    
    # define the augmentation of covariates,
    # which is the combination of violation space and baseline covariates
    getDesign_rf <- function(X, knots=5) {
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
    Cov.aug = cbind(vio.space, getDesign_rf(X))

    # Treatment model fitting
    forest <- TSCI.RF.fit(D,Z,X,
                          num.trees = 500,
                          mtry = seq(round(P/3),round(2*P/3),by = 1),
                          max.depth = 0,
                          min.node.size = c(5,10,20),
                          split.prop = 2/3)
    A1.ind <- forest$A1.ind
    # Compute weight matrix
    weight <- TSCI.RF.weight(forest$nodes.A1)

    # Selection
    res_tsci <- TSCI.RF.Selection(Y,D,Cov.aug,
                                  A1.ind,
                                  weight = weight,
                                  Q = Q,
                                  rm.ind = rm.ind,
                                  intercept = T,
                                  layer = T,
                                  str.thol = 10,
                                  alpha = 0.05)

    # TSCI-basis --------------------------------------------------------------

    # model fitting
    basis.fit <- TSCI.basis.fit(D, cbind(Z,X))
    D.rep <- basis.fit$D.rep
    M <- basis.fit$M
    knot <- basis.fit$knot
    
    # basis for X
    X_basis = getDesign_rf(X)
    rm_cols = seq(1, ncol(X_basis), by = 7)
    X_basis = X_basis[,-rm_cols]
    W = cbind(Z, X_basis)
    
    res_tsciba_hetero <- TSCI.basis.Selection.hetero(Y, D, W, 
                                                     D.rep, 
                                                     knot, 
                                                     M, 
                                                     Q, 
                                                     vio.space, 
                                                     intercept = T, 
                                                     str.thol = 10)

    # TSLS --------------------------------------------------------------------
    
    res_tsls = ivmodel(Y, D, Z, X, heteroSE = T)
    coef_tsls = coef(res_tsls)["TSLS", "Estimate"]
    sd_tsls = coef(res_tsls)["TSLS", "Std. Error"]
    # concentration parameter
    stage1 = lm(D ~ Z + X)
    delta_hat = resid(stage1)
    nume = sum(resid(lm(predict(stage1) ~ X))^2)
    concen = nume / mean(delta_hat^2)

    # DML ---------------------------------------------------------------------

    ml_l = lrn("regr.ranger", num.trees = 500, max.depth = 0)
    ml_m = ml_l$clone()
    ml_r = ml_l$clone()
    dt = data.table(cbind(Y, D, Z, X))
    colnames(dt) = c("y","d","z",paste0("x",1:p))
    obj_dml_data = DoubleMLData$new(dt, y_col = "y", d_cols = "d", 
                                    x_cols = paste0("x",1:p), z_cols = "z")
    dml_pliv_obj = DoubleMLPLIV$new(obj_dml_data, ml_l, ml_m, ml_r)
    P = p + 1
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
    
    # return ------------------------------------------------------------------
    
    returnValue = list(res_tsci = res_tsci
                       , res_tsciba_hetero = res_tsciba_hetero
                       , res_tsls = c(coef_tsls, sd_tsls, concen)
                       , res_dml = c(beta_DMLorg, sd_DMLorg)
    )
    returnValue
}

nsim = 20
total = NULL
for (i in 1:nsim) {
    print(i)
    result = sim_func(a)
    total = rbind(total, result)
}
filename = paste0("results-inter-strongHidden/results"
                  , "-a", a
                  , "-round", round, ".rds")
saveRDS(total, filename)

