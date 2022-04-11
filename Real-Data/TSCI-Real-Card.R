# card data
library(ivmodel)
data(card.data)
source("Source-RF-hetero2.R",encoding="UTF-8")


nsim = 10
results.list = rep(list(NA),nsim)
results.list2 = rep(list(NA),nsim)
Xname1=c("exper", "expersq", "black", "south", "smsa")
Xname2=c("exper", "expersq", "black", "south", "smsa", "reg661","reg662", "reg663", "reg664", "reg665", "reg666", "reg667","reg668", "smsa66")

for (i in 1:20) {
  print(i)
  ind = sample(1:3010,3010,replace = FALSE)
  Y=as.matrix(card.data[ind,"lwage"])
  D=as.matrix(card.data[ind,"educ"])
  Z=as.matrix(card.data[ind,"nearc4"])
  
  # X1=card.data[,Xname1]
  # foo1 = ivmodel(Y=Y,D=D,Z=Z,X=X1)
  # summary(foo1)
  # 
  # rf1 = TSCI.RF(Y,D,Z,X1,vio.space = Z,mtry = 1:2,min.node.size = c(5,10,20,30))
  
  
  X2=as.matrix(card.data[ind,Xname2])
  # foo1 = ivmodel(Y=Y,D=D,Z=Z,X=X2)
  # summary(foo1)
  
  inter = matrix(NA,nrow(Z),5)
  for (j in 1:5) {
    inter[,j] = Z*X2[,j]
  }
  vio.space = list(inter,Z)
  vio.space2 = list(cbind(inter,Z))
  rf.card = TSCI.RF(Y,D,Z,X2,vio.space = vio.space,mtry = 1:(1+length(Xname2)),min.node.size = c(5,10,20,30),num.trees = 500)
  rf.card2 = TSCI.RF(Y,D,Z,X2,vio.space = vio.space2,mtry = 1:(1+length(Xname2)),min.node.size = c(5,10,20,30),num.trees = 500)
  
  results.list[[i]] = rf.card
  results.list2[[i]] = rf.card2
}




# # full RF
# source("Source-otherRF-homo.R",encoding="UTF-8")
# p = ncol(X2)+1
# Y = as.matrix(Y)
# X2 = as.matrix(X2)
# D = as.matrix(D)
# Z = as.matrix(Z)
# 
# mtry = seq(round((p+1)/3), round(2*(p+1)/3), by=1)
# ### set max.depth and min.node.size for tuning
# ### larger max.depth and smaller min.node.size means more complex trees
# max.depth <- 0; min.node.size <- c(5,10,20)
# 
# 
# forest.cov = cbind(Z,X2)
# forest.cov = as.matrix(forest.cov)
# 
# forest.full <- TSRF.full(forest.cov,D,mtry=mtry,max.depth=max.depth,min.node.size=min.node.size,num.trees = 500)
# weight.full <- TSCI.RF.weight(forest.full$nodes)
# full.outputs <- TSRF.stat.full(Y,D,forest.cov,weight.full)
# 
# full.outputs2 <- TSRF.stat.full(Y,D,X2,weight.full)
# 
# full.outputs$betaHat
# full.outputs$sd
# full.outputs$iv.str/full.outputs$SigmaSqD


