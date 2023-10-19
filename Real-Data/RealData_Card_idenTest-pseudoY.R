library(ivmodel)
data(card.data)
source("Source-RF-hetero2.R", encoding="UTF-8")

# parameter change --------------------------------------------------------

round = 1  # {1, 2, ..., 25}

# run invalidity test -----------------------------------------------------

nsplits = 20
Xname1=c("exper", "expersq", "black", "south", "smsa", "smsa66")
Xname2=c("exper", "expersq", "black", "south", "smsa", "reg661", "reg662"
         , "reg663", "reg664", "reg665", "reg666", "reg667","reg668", "smsa66")

total = NULL
for (i in 1:nsplits) {
    print(i)
    
    ### run TSCI
    ind = sample(1:nrow(card.data), nrow(card.data), replace=FALSE)
    Y = as.matrix(card.data[ind,"lwage"])
    D = as.matrix(card.data[ind,"educ"])
    Z = as.matrix(card.data[ind,"nearc4"])
    X1 = as.matrix(card.data[ind,Xname1])
    X2 = as.matrix(card.data[ind,Xname2])
    # build violation space
    V1 = matrix(NA, nrow(Z), ncol(X1)+1)
    V1[,1] = Z
    for (j in 1:ncol(X1)) {
        V1[,j+1] = Z * X1[,j]
    }
    V2 = matrix(NA, nrow(Z), ncol(X2)-ncol(X1))
    for (j in 1:ncol(V2)) {
        V2[,j] = Z * X2[,(5+j)]
    }
    vio_space = list(V2, V1)
    # TSCI
    rf_card = TSCI.RF(Y,D,Z,X2
                      ,vio.space=vio_space
                      ,mtry = 1:(1+length(Xname2))
                      ,min.node.size = c(5,10,20,30)
                      ,num.trees = 500)
    
    ### test nonlinearity
    betaHat = rf_card$Coef.robust[2]
    qcomp = rf_card$q.comp
    # split samples
    ind_test = sample(1:nrow(card.data), round(nrow(card.data)*2/3), replace=FALSE)
    Y_train = as.matrix(card.data[-ind_test,"lwage"])
    D_train = as.matrix(card.data[-ind_test,"educ"])
    Z_train = as.matrix(card.data[-ind_test,"nearc4"])
    X1_train = as.matrix(card.data[-ind_test,Xname1])
    X_train = as.matrix(card.data[-ind_test,Xname2])
    Y_test = as.matrix(card.data[ind_test,"lwage"])
    D_test = as.matrix(card.data[ind_test,"educ"])
    Z_test = as.matrix(card.data[ind_test,"nearc4"])
    X1_test = as.matrix(card.data[ind_test,Xname1])
    X_test = as.matrix(card.data[ind_test,Xname2])
    V1_train = matrix(NA, nrow(Z_train), ncol(X1)+1)
    V1_train[,1] = Z_train
    for (j in 1:ncol(X1)) {
        V1_train[,j+1] = Z_train * X1_train[,j]
    }
    V1_test = matrix(NA, nrow(Z_test), ncol(X1)+1)
    V1_test[,1] = Z_test
    for (j in 1:ncol(X1)) {
        V1_test[,j+1] = Z_test * X1_test[,j]
    }
    V2_train = matrix(NA, nrow(Z_train), ncol(X2)-ncol(X1))
    for (j in 1:ncol(V2)) {
        V2_train[,j] = Z_train * X_train[,(5+j)]
    }
    V2_test = matrix(NA, nrow(Z_test), ncol(X2)-ncol(X1))
    for (j in 1:ncol(V2)) {
        V2_test[,j] = Z_test * X_test[,(5+j)]
    }
    psuedo_Y_train = Y_train - betaHat * D_train
    psuedo_Y_test = Y_test - betaHat * D_test
    # regression with V_0
    reg = lm(psuedo_Y_train ~ X_train)
    coef1 = coef(reg)
    pred1 = cbind(1, X_test) %*% coef1
    resid_reg = psuedo_Y_test - pred1
    sigmaSq1 = mean(resid_reg^2)
    # regression with V_1
    reg = lm(psuedo_Y_train ~ X_train + V1_train)
    coef2 = coef(reg)
    pred2 = cbind(1, X_test, V1_test) %*% coef2
    resid_reg = psuedo_Y_test - pred2
    sigmaSq2 = mean(resid_reg^2)
    # regression with V_2
    reg = lm(psuedo_Y_train ~ X_train + V1_train + V2_train)
    coef3 = coef(reg)
    pred3 = cbind(1, X_test, V1_test, V2_test) %*% coef3
    resid_reg = psuedo_Y_test - pred3
    sigmaSq3 = mean(resid_reg^2)
    # random forests with any fitted V
    ranger_Data_train = cbind(psuedo_Y_train, Z_train, X_train)
    ranger_Data_test = cbind(psuedo_Y_test, Z_test, X_test)
    colnames(ranger_Data_train) = c("psuedo_Y_train", paste0("W", 1:(1+ncol(X_train))))
    colnames(ranger_Data_test) = c("psuedo_Y_test", paste0("W", 1:(1+ncol(X_test))))
    num.trees = 500
    mtry = 1:(1+ncol(X_train))
    max.depth = 0
    min.node.size = c(5,10,20,30)
    params.grid = expand.grid(
        num.trees = num.trees,
        mtry = mtry,
        max.depth = max.depth,
        min.node.size = min.node.size
    )
    sigmaSq_rf = 10000
    for (i in 1:nrow(params.grid)) {
        temp_forest = ranger(psuedo_Y_train~., data = ranger_Data_train,
                             num.trees=params.grid$num.trees[i],
                             mtry=params.grid$mtry[i],
                             max.depth = params.grid$max.depth[i],
                             min.node.size = params.grid$min.node.size[i],
                             importance = "impurity"
        )
        temp_pred = predict(temp_forest, data=ranger_Data_test)$predictions
        resid_temp = psuedo_Y_test - temp_pred
        sigmaSq_rf_temp = mean(resid_temp^2)
        if (sigmaSq_rf_temp <= sigmaSq_rf) {
            pred_rf = temp_pred
            sigmaSq_rf = sigmaSq_rf_temp
            param = params.grid[i,]
        }
    }
    # record results
    results = list(coef1 = coef1
                   , coef2 = coef2
                   , coef3 = coef3
                   , pred = cbind(pred1, pred2, pred3, pred_rf, psuedo_Y_test)
                   , sigmaSq = c(sigmaSq1, sigmaSq2, sigmaSq3, sigmaSq_rf)
                   , forest.param = param)
    total = rbind(total, results)
}
filename = paste0("results-RealData_Card_invalidityTest-pseudoY-round", round, ".rds")
saveRDS(total, filename)



