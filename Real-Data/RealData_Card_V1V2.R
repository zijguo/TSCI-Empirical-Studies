library(ivmodel)
data(card.data)
source("Source-RF-hetero2.R", encoding="UTF-8")

# parameter change --------------------------------------------------------

round = 1  # {1, 2, ..., 25}

# run algorithm -----------------------------------------------------------

nsplits = 20
results.list = NULL
Xname1=c("exper", "expersq", "black", "south", "smsa", "smsa66")
Xname2=c("exper", "expersq", "black", "south", "smsa", "reg661", "reg662"
         , "reg663", "reg664", "reg665", "reg666", "reg667","reg668", "smsa66")

for (i in 1:nsplits) {
    print(i)
    ind = sample(1:nrow(card.data), nrow(card.data), replace=FALSE)
    Y = as.matrix(card.data[ind,"lwage"])
    D = as.matrix(card.data[ind,"educ"])
    Z = as.matrix(card.data[ind,"nearc4"])
    X1 = as.matrix(card.data[ind,Xname1])
    X2 = as.matrix(card.data[ind,Xname2])
    
    ### build violation spaces
    V1 = matrix(NA, nrow(Z), ncol(X1)+1)
    V1[,1] = Z
    for (j in 1:ncol(X1)) {
        V1[,j+1] = Z * X1[,j]
    }

    # first V2 (include other Z * X)
    V2_1 = matrix(NA, nrow(Z), ncol(X2)-ncol(X1))
    for (j in 1:ncol(V2_1)) {
        V2_1[,j] = Z * X2[,(5+j)]
    }
    
    # second V2 (only include Z * exper^3)
    V2_2 = Z * X2[,1]^3
    
    # third V2 (include Z * two-way interaction of 6 most important X, excluding exper^3)
    V2_3 = NULL
    for (k in 1:(ncol(X1)-1)) {
        for (l in (k+1):ncol(X1)) {
            if (Xname1[k]=='exper' & Xname1[l]=='expersq') next
            V2_3 = cbind(V2_3, Z * X1[,k] * X1[,l])
        }
    }
    V2_3 = V2_3[, !duplicated(t(cbind(V1, V2_3)))[(ncol(V1)+1):ncol(cbind(V1, V2_3))]]
    
    # fourth V2 (include Z * two-way interaction of 6 most important X - including exper^3)
    V2_4 = NULL
    for (k in 1:(ncol(X1)-1)) {
        for (l in (k+1):ncol(X1)) {
            # if (Xname1[k]=='exper' & Xname1[l]=='expersq') next
            V2_4 = cbind(V2_4, Z * X1[,k] * X1[,l])
        }
    }
    V2_4 = V2_4[, !duplicated(t(cbind(V1, V2_4)))[(ncol(V1)+1):ncol(cbind(V1, V2_4))]]
    
    vio_space_1 = list(V2_1, V1)
    vio_space_2 = list(V2_2, V1)
    vio_space_3 = list(V2_3, V1)
    vio_space_4 = list(V2_4, V1)
    
    ### run TSCI
    rf_card_1 = TSCI.RF(Y,D,Z,X2
                         ,vio.space=vio_space_1
                         ,mtry = 1:(1+length(Xname2))
                         ,min.node.size = c(5,10,20,30)
                         ,num.trees = 500)
    rf_card_2 = TSCI.RF(Y,D,Z,X2
                          ,vio.space=vio_space_2
                          ,mtry = 1:(1+length(Xname2))
                          ,min.node.size = c(5,10,20,30)
                          ,num.trees = 500)
    rf_card_3 = TSCI.RF(Y,D,Z,X2
                        ,vio.space=vio_space_3
                        ,mtry = 1:(1+length(Xname2))
                        ,min.node.size = c(5,10,20, 30)
                        ,num.trees = 500)
    rf_card_4 = TSCI.RF(Y,D,Z,X2
                        ,vio.space=vio_space_4
                        ,mtry = 1:(1+length(Xname2))
                        ,min.node.size = c(5,10,20,30)
                        ,num.trees = 500)
    result = list(rf_card_1, rf_card_2, rf_card_3, rf_card_4)
    results.list = rbind(results.list, result)
}

filename = paste0('results-RealData_Card_largerVio-round', round, '.rds')
saveRDS(results.list, filename)

