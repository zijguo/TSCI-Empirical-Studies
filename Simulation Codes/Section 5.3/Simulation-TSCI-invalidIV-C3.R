library(ivmodel)
library(fda)
library(randomForest)
library(DoubleML)
library(mlr3)
library(mlr3learners)
library(mlr3tuning)
library(data.table)
source("Source-Basis-hetero.R")
source("Source-RF-hetero2.R")
source("Source-DML_helpers.R")

# parameter change -------------------------------------------------------

round = 1  # {1, 2, ..., 25}
# nonliearity degree
a = 0  # {0, 0.1, 0.2, ..., 1.2}

# simulation --------------------------------------------------------------

sim_func <- function(a) {
    # generate data -----------------------------------------------------------
    
    # dimension
    p = 5
    n = 3000
    # X and Z
    mu_Xstar = rep(0, p + 1)
    cov_Xstar = stats::toeplitz(0.5^c(0:p))
    Xstar = MASS::mvrnorm(n, mu_Xstar, cov_Xstar)
    X = pnorm(Xstar)[,-1]
    Z = (pnorm(Xstar)[,1] - 0.5) * 4
    # hidden confounders
    h <- 2*sin(X[,1]) + rnorm(n, 0, 1)
    # D
    f_func = function(Z, X, a) {
        Z + a*(sin(2*pi*Z) + 1.5*cos(2*pi*Z)) - 0.3*rowSums(X)
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
    res_tsci <- TSCI.RF.Selection(Y,D,Cov.aug,A1.ind,weight=weight,rm.ind=rm.ind,Q=Q,intercept=T,layer=T,str.thol=10,alpha=0.05)

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
    concen_tsls = nume / mean(delta_hat^2)

    # DML (use detailed steps) ------------------------------------------------
    
    ### preprocess data
    aa <- as.matrix(Z)
    xx <- as.matrix(D)
    yy <- as.matrix(Y)
    ww <- as.matrix(X)
    colnames(aa) = "aa"
    colnames(xx) = "xx"
    colnames(yy) = "yy"
    # colnames(ww) = "ww"
    colnames(ww) = paste0("ww", 1:ncol(ww))
    d <- ncol(xx)
    xx_colnames <- apply(rbind(rep("b", d), seq_len(d)), 2,
                         function(x) paste(x[1], x[2], sep = ""))
    
    ### specify ML methods
    P = length(colnames(ww)) + 1
    cond_func_all = get_condexp_funcs(
        cond_method = c("forest", "forest", "forest"),
        params = list(list(min.node.size = c(5,10,20), mtry = seq(round(P/3), round(2*P/3), by=1)),
                      list(min.node.size = c(5,10,20), mtry = seq(round(P/3), round(2*P/3), by=1)),
                      list(min.node.size = c(5,10,20), mtry = seq(round(P/3), round(2*P/3), by=1)))
    )
    
    ### fit the ML models (include cross-fitting)
    K = 2
    n_reorder <- sample(seq_len(n), n, replace = FALSE)
    folds = cut(n_reorder, breaks = K, labels = FALSE)
    all_residuals = all_residuals_gW = vector(mode = "list", length = K)
    for (kk in seq_len(K)) {
        I <- seq_len(n)[folds == kk] # test set: evaluate conditional expectations
        Ic <- setdiff(seq_len(n), I) # training set: estimate conditional expectations
        
        # train conditional expectations with Ic, and evaluate residuals on I.
        # cond_func_all = learning method for conditional expectations
        # params = additional parameters to estimate conditional expectations
        residuals_samplesplit_hat <-
            residuals_samplesplit(aa = aa, ww = ww,
                                  xx = xx, yy = yy,
                                  I = I, Ic = Ic,
                                  cond_func_all = cond_func_all,
                                  params = NULL)
        
        all_residuals[[kk]] <- residuals_samplesplit_hat[c("rA", "rX", "rY")]
        all_residuals_gW[[kk]] <- residuals_samplesplit_hat[c("rX_gW", "rY_gW")]
    } # end for (kk in seq_len(K))
    
    ### generate the point estimator
    # there are 2 decomposition ways for the estimator
    # the first one is beta = gamma*rA_sq*Gamma / gamma*rA_sq*gamma (we can see the weights of 1,2,..,K) -> to compare with TSCI
    # the second one is beta = Gamma_nume / gamma_nume                                                   -> to compare between TSLS & DML
    K <- length(all_residuals_gW)
    dim_rX <- dim(all_residuals_gW[[1]]$rX_gW)
    d <- dim_rX[2]
    n <- dim_rX[1]
    mat <- matrix(0, nrow = d, ncol = d)
    vec <- matrix(0, nrow = d, ncol = 1)
    beta_nume_DML = beta_deno_DML = rep(0, K)
    gamma_DML = Gamma_DML = rA_sq_DML = rep(0, K)
    rX_concen = rA_concen = NULL
    for (k in seq_len(K)) {
        mat <- mat + crossprod(all_residuals_gW[[k]]$rX_gW, all_residuals_gW[[k]]$rX_gW) / n
        vec <- vec + crossprod(all_residuals_gW[[k]]$rX_gW, all_residuals_gW[[k]]$rY_gW) / n
        rA = all_residuals[[k]]$rA
        rX = all_residuals[[k]]$rX
        rY = all_residuals[[k]]$rY
        beta_deno_DML[k] = sum(rA*rX)
        beta_nume_DML[k] = sum(rA*rY)
        gamma_DML[k] = qr.solve(rA, rX)
        Gamma_DML[k] = qr.solve(rA, rY)
        rA_sq_DML[k] = sum(rA^2)
        rX_concen = c(rX_concen, rX)
        rA_concen = c(rA_concen, rA)
    }
    beta_DML = qr.solve(mat, vec)
    
    ### generate the variance estimator
    var_DML = sigma2_DML(all_residuals = all_residuals,
                         betahat = beta_DML)
    sd_DML = sqrt(var_DML)
    
    ### concentration parameter (residuals)
    stage1 = lm(rX_concen~rA_concen-1)
    nume = sum(predict(stage1)^2)
    concen_DML = nume / mean(resid(stage1)^2)
   
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
    
    returnValue = list(res_tsci = res_tsci,
                       res_tsciba_hetero = res_tsciba_hetero,
                       res_tsls = c(coef_tsls, sd_tsls, concen_tsls),
                       res_dml = c(beta_DMLorg, sd_DMLorg, concen_DML)
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
filename = paste0("results-sin-littleHidden/results"
                  , "-a", a
                  , "-round", round, ".rds")
saveRDS(total, filename)

