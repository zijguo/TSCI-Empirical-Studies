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
# IV strength
a = 1  # {0.6, 0.8, ..., 2}

# simulation --------------------------------------------------------------

sim_func = function(a) {
    ## Generate some data
    # true linear parameter
    beta0 <- 1
    n <- 3000
    # observed confounder
    X <- pi * runif(n, -1, 1)
    # instrument
    # Z <- 3 * tanh(2 * X - 1) + rnorm(n, 0, 1)
    Z = 0.6 * X + rnorm(n, 0, 1)
    # unobserved confounder
    h <- 2 * sin(X) + rnorm(n, 0, 1)
    # linear covariate
    # D <- -a * abs(Z) - h - 2 * tanh(X) + rnorm(n, 0, 1)
    D = a*Z^2/2 - h - 2 * tanh(X) + rnorm(n, 0, 1)
    # response
    Y <- beta0 * D - 3 * cos(pi * 0.25 * h) + 0.5 * X ^ 2 + rnorm(n, 0, 1)

    # DML ---------------------------------------------------------------------

    ml_l = lrn("regr.ranger", num.trees = 500, max.depth = 0)
    ml_m = ml_l$clone()
    ml_r = ml_l$clone()
    dt = data.table(cbind(Y, D, Z, X))
    colnames(dt) = c("y","d","z","x")
    obj_dml_data = DoubleMLData$new(dt, y_col = "y", d_cols = "d", 
                                    x_cols = "x", z_cols = "z")
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
    
    # DML (use detailed steps) ------------------------------------------------
    
    # ### preprocess data
    # aa <- as.matrix(Z)
    # xx <- as.matrix(D)
    # yy <- as.matrix(Y)
    # ww <- as.matrix(X)
    # colnames(aa) = "aa"
    # colnames(xx) = "xx"
    # colnames(yy) = "yy"
    # # colnames(ww) = "ww"
    # colnames(ww) = paste0("ww", 1:ncol(ww))
    # d <- ncol(xx)
    # xx_colnames <- apply(rbind(rep("b", d), seq_len(d)), 2,
    #                      function(x) paste(x[1], x[2], sep = ""))
    # 
    # ### specify ML methods
    # P = length(colnames(ww)) + 1
    # cond_func_all = get_condexp_funcs(
    #     cond_method = c("forest", "forest", "forest"),
    #     params = list(list(min.node.size = c(5,10,20), mtry = seq(round(P/3), round(2*P/3), by=1)),
    #                   list(min.node.size = c(5,10,20), mtry = seq(round(P/3), round(2*P/3), by=1)),
    #                   list(min.node.size = c(5,10,20), mtry = seq(round(P/3), round(2*P/3), by=1)))
    # )
    # 
    # ### fit the ML models (include cross-fitting)
    # K = 2
    # n_reorder <- sample(seq_len(n), n, replace = FALSE)
    # folds = cut(n_reorder, breaks = K, labels = FALSE)
    # all_residuals = all_residuals_gW = vector(mode = "list", length = K)
    # for (kk in seq_len(K)) {
    #     I <- seq_len(n)[folds == kk] # test set: evaluate conditional expectations
    #     Ic <- setdiff(seq_len(n), I) # training set: estimate conditional expectations
    #     
    #     # train conditional expectations with Ic, and evaluate residuals on I.
    #     # cond_func_all = learning method for conditional expectations
    #     # params = additional parameters to estimate conditional expectations
    #     residuals_samplesplit_hat <-
    #         residuals_samplesplit(aa = aa, ww = ww,
    #                               xx = xx, yy = yy,
    #                               I = I, Ic = Ic,
    #                               cond_func_all = cond_func_all,
    #                               params = NULL)
    #     
    #     all_residuals[[kk]] <- residuals_samplesplit_hat[c("rA", "rX", "rY")]
    #     all_residuals_gW[[kk]] <- residuals_samplesplit_hat[c("rX_gW", "rY_gW")]
    # } # end for (kk in seq_len(K))
    # 
    # ### generate the point estimator
    # # there are 2 decomposition ways for the estimator
    # # the first one is beta = gamma*rA_sq*Gamma / gamma*rA_sq*gamma (we can see the weights of 1,2,..,K) -> to compare with TSCI
    # # the second one is beta = Gamma_nume / gamma_nume                                                   -> to compare between TSLS & DML
    # K <- length(all_residuals_gW)
    # dim_rX <- dim(all_residuals_gW[[1]]$rX_gW)
    # d <- dim_rX[2]
    # n <- dim_rX[1]
    # mat <- matrix(0, nrow = d, ncol = d)
    # vec <- matrix(0, nrow = d, ncol = 1)
    # beta_nume_DML = beta_deno_DML = rep(0, K)
    # gamma_DML = Gamma_DML = rA_sq_DML = rep(0, K)
    # rX_concen = rA_concen = NULL
    # for (k in seq_len(K)) {
    #     mat <- mat + crossprod(all_residuals_gW[[k]]$rX_gW, all_residuals_gW[[k]]$rX_gW) / n
    #     vec <- vec + crossprod(all_residuals_gW[[k]]$rX_gW, all_residuals_gW[[k]]$rY_gW) / n
    #     rA = all_residuals[[k]]$rA
    #     rX = all_residuals[[k]]$rX
    #     rY = all_residuals[[k]]$rY
    #     beta_deno_DML[k] = sum(rA*rX)
    #     beta_nume_DML[k] = sum(rA*rY)
    #     gamma_DML[k] = qr.solve(rA, rX)
    #     Gamma_DML[k] = qr.solve(rA, rY)
    #     rA_sq_DML[k] = sum(rA^2)
    #     rX_concen = c(rX_concen, rX)
    #     rA_concen = c(rA_concen, rA)
    # }
    # beta_DML = qr.solve(mat, vec)
    # 
    # ### generate the variance estimator
    # var_DML = sigma2_DML(all_residuals = all_residuals,
    #                      betahat = beta_DML)
    # sd_DML = sqrt(var_DML)
    # 
    # ### concentration parameter (residuals)
    # stage1 = lm(rX_concen~rA_concen-1)
    # nume = sum(predict(stage1)^2)
    # concen_DML = nume / mean(resid(stage1)^2)
    
    # regsDML -----------------------------------------------------------------
    
    # res_regsdml = regsdml(a = Z, w = X, x = D, y = Y,
    #                       do_regsDML = T,
    #                       # gamma = exp(seq(-4, 1, length.out = 4)),
    #                       do_regDML_all_gamma = T,
    #                       S = 1, K = 2,
    #                       cond_method = c("forest","forest","forest"),
    #                       params = list(list(min.node.size = c(5,10,20), mtry = seq(round(P/3), round(2*P/3), by=1)),
    #                                     list(min.node.size = c(5,10,20), mtry = seq(round(P/3), round(2*P/3), by=1)),
    #                                     list(min.node.size = c(5,10,20), mtry = seq(round(P/3), round(2*P/3), by=1))))
    # beta_regsdml = res_regsdml$regsDML_statistics$beta_regsDML
    # sd_regsdml = res_regsdml$regsDML_statistics$sd_regsDML
    
    # TSLS --------------------------------------------------------------------
    
    # same as the above, we use 2 decomposition ways
    # beta_tsls = coef(ivreg(Y ~ D+X | Z+X))[2]
    model_tsls = ivmodel(Y, D, Z, X, heteroSE = T)
    beta_tsls = coef(model_tsls)["TSLS", "Estimate"]
    
    ### variance estimator
    eps_hat = resid(lm(Y-D*beta_tsls~X))
    # sd_tsls = sqrt(mean(eps_hat^2) * rA_sq_tsls / beta_deno_tsls^2)
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
        # DML
        # beta_DML = beta_DML,
        # gamma_DML = gamma_DML,
        # Gamma_DML = Gamma_DML,
        # rA_sq_DML = rA_sq_DML,
        # beta_nume_DML = beta_nume_DML,
        # beta_deno_DML = beta_deno_DML,
        # sd_DML = sd_DML,
        # concen_DML = concen_DML,
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

total = NULL
nsim = 20
for (i in 1:nsim) {
    print(i)
    result = sim_func(a)
    total = rbind(total, result)
}
filename = paste0("results-Zsq/results",
                  "-a", a,
                  "-round", round, ".rds")
saveRDS(total, filename)

