# simulation of other methods in homoskedastic settings


library(MASS)
source('~/Source-otherRF-homo.R', encoding = 'UTF-8')
source("~/Source-Basis-homo.R",encoding = "UTF-8")


# function to generate covariance matrix
A1gen=function(rho,p){
  A1=matrix(0,p,p)
  for(i in 1:p){
    for(j in 1:p){
      A1[i,j]=rho^(abs(i-j))
    }
  }
  A1
}



###### dimension change across 20
p = 20
####please set n = 5000
n = 1000
### setting, change across 1,2
f.index = 1
##### change the interaction across 0, 0.5, 1
inter.val = 1
#### violation index, change across 1,2
vio.index = 2
### change round from 1 to 20
round = 1
#### a denotes the IV strength, set as 1
a = 1
##### tau denotes the violation strength, set as 1
tau = 1
#### the number of simulations, make sure there is 500 in total
#### if you specify round to be 1,2,...,10, then set nsim as 50
nsim = 10


######### DO NOT MODIFY BELOW #########
f_1 = function(x){x+a*(x^2+0.125*x^4) -25/12}
# f_2 = function(x){exp(2*a+x+0.5*x^3)+a*(x^2+0.5*x^4)}
f_2 = function(x){a*(1*sin(2*pi*x) + 1.5*cos(2*pi*x))}
rho1=0.5
Cov=(A1gen(rho1,p+1))
mu=rep(0,p+1)
beta=1
alpha=as.matrix(rep(-0.3,p))
gamma=as.matrix(rep(0.2,p))
inter=as.matrix(c(rep(inter.val,5),rep(0,p-5)))

# number of violation spaces
Q = 4
# estimator names
estimator.names = paste("q",0:(Q-1),sep="")
### the fixed violation space for RF
Coef.matrix.basis=matrix(NA,nrow=nsim,ncol=Q)
sd.matrix.basis=matrix(NA,nrow=nsim,ncol=Q)
colnames(Coef.matrix.basis) = colnames(sd.matrix.basis) = estimator.names


# IV strength test
iv.str = iv.thol = matrix(NA,nrow=nsim,ncol=Q)
colnames(iv.str) = colnames(iv.thol) = estimator.names


# # variance of epsilon
SigmaSqY = matrix(NA, nsim, Q)
colnames(SigmaSqY) = estimator.names
# variance of delta
SigmaSqD = matrix(NA, nsim, 1)
SigmaSqY.Qmax = SigmaSqD


### selection
Qmax  = matrix(NA,nsim,1)
q.comp = q.robust = SigmaSqY.Qmax = Qmax
Coef.robust = matrix(NA,nsim,2)
sd.robust = matrix(NA,nsim,2)
colnames(Coef.robust) = colnames(sd.robust) = c("Basis-comp","Basis-robust")

### oracle estimator
Coef.oracle = sd.oracle = matrix(NA,nsim,1)
colnames(Coef.oracle) = colnames(sd.oracle) = c("Basis")


### naive/plain estimators
Coef.naive <- sd.naive <- matrix(NA,nsim,2)
colnames(Coef.naive) <- colnames(sd.naive) <- c("naiveRF","RF-full")


### TSLS estimator
Coef.TSLS <- sd.TSLS <- matrix(NA,nsim,1)

### weak IV problem
run.OLS = weak.iv = matrix(NA,nsim,1)
colnames(run.OLS) = colnames(weak.iv) = "Basis"


for(i in 1:nsim) {
  print(i)
  #### generate the data
  mu.error=rep(0,2)
  Cov.error=matrix(c(1,0.5,0.5,1),2,2)
  Error=mvrnorm(n, mu.error, Cov.error)
  W.original=mvrnorm(n, mu, Cov)
  W=pnorm(W.original)
  Z=W[,1]
  Z = 4*(Z-0.5)
  X=W[,-1]
  ###### generate the data for the treatment variable D
  if(f.index==1){
    D=f_1(Z)+X%*%alpha+Z*X%*%inter+Error[,1]
  }
  if(f.index==2){
    D=f_2(Z)+X%*%alpha+Z*X%*%inter+Error[,1]
  }
  ####### generate the outcome variable
  # if(vio.index==0){
  #   Y=D*beta+ X%*% gamma+Error[,2]
  # }
  if(vio.index==1){
    Y=D*beta+ tau*Z+ X%*%gamma+Error[,2]
  }
  if(vio.index==2){
    Y=D*beta+ tau*(Z^2+Z-1)+ X%*% gamma+Error[,2] # difficult if no Z
  }

  # ### TSLS Estimators
  # TSLS <- ivreg(Y~D+Z+X|Z+X)
  # Coef.TSLS[i,1] <- coef(TSLS)[2]
  # sd.TSLS[i,1] <- summary(TSLS)$coefficients[2,2]
  # 
  # # basis appraoch
  # outputs.basis = TSCI.basis(Y,D,Z,X)
  # if (length(outputs.basis$Coef.vec)>=Q) {
  #   Coef.matrix.basis[i,] <- outputs.basis$Coef.vec
  #   sd.matrix.basis[i,] <- outputs.basis$sd.vec
  #   SigmaSqY[i,] <- outputs.basis$SigmaSqY
  #   iv.str[i,] = outputs.basis$iv.str
  #   iv.thol[i,] = outputs.basis$iv.thol
  #   
  # } else {
  #   Coef.matrix.basis[i,] <- c(outputs.basis$Coef.vec, rep(NA,Q-length(outputs.basis$Coef.vec)))
  #   sd.matrix.basis[i,] <- c(outputs.basis$sd.vec, rep(NA,Q-length(outputs.basis$sd.vec)))
  #   SigmaSqY[i,] <- c(outputs.basis$SigmaSqY,rep(NA,Q-length(outputs.basis$SigmaSqY)))
  #   iv.str[i,] = c(outputs.basis$iv.str,rep(NA,Q-length(outputs.basis$iv.str)))
  #   iv.thol[i,] = c(outputs.basis$iv.thol,rep(NA,Q-length(outputs.basis$iv.thol)))
  # }
  # Coef.robust[i,] <- outputs.basis$Coef.robust
  # sd.robust[i,] <- outputs.basis$sd.robust
  # Coef.oracle[i,1] <- outputs.basis$Coef.vec[vio.index+1]
  # sd.oracle[i,1] <- outputs.basis$sd.vec[vio.index+1]
  # 
  # SigmaSqY.Qmax[i,] <- outputs.basis$SigmaSqY.Qmax
  # SigmaSqY[i,] <- outputs.basis$SigmaSqY
  # SigmaSqD[i,] <- outputs.basis$SigmaSqD
  # 
  # 
  # q.comp[i,] = outputs.basis$q.comp
  # q.robust[i,] = outputs.basis$q.robust
  # Qmax[i,] = outputs.basis$Q.max
  # run.OLS[i,] = outputs.basis$run.OLS
  # weak.iv[i,] = outputs.basis$weak.iv
  # 
  # 
  # 
  # # other RF methods
  mtry = seq(round((p+1)/3), round(2*(p+1)/3), by=1)
  ### set max.depth and min.node.size for tuning
  ### larger max.depth and smaller min.node.size means more complex trees
  max.depth <- 0; min.node.size <- c(5,10,20)
  # 
  # forest.2 <- TSCI.RF.fit(D,Z,X,mtry=mtry,max.depth=max.depth,min.node.size=min.node.size,split.prop=2/3,num.trees = 200)
  # A1.ind <- forest.2$A1.ind
  # # weight matrix
  # weight.2 <- TSCI.RF.weight(forest.2$nodes.A1)
  # rm(forest.2)
  # 
  Cov.aug <- W[,-1]
  for (q in 1:(Q-1)) {
    Cov.aug<-cbind(Z^q,Cov.aug)
  }

  # 
  # 
  # naiveRF.outputs <- naiveRF.stat(Y, D, Cov.aug[,-(1:(Q-1-vio.index))], A1.ind, weight.2)
  # Coef.naive[i,1] <- naiveRF.outputs$betaHat
  # sd.naive[i,1] <- naiveRF.outputs$sd
  # rm(weight.2)
  
  
  ### random forest with full data
  forest.cov = cbind(Z,X)
  forest.full <- TSRF.full(forest.cov,D,mtry=mtry,max.depth=max.depth,min.node.size=min.node.size)
  weight.full <- TSCI.RF.weight(forest.full$nodes)
  full.outputs <- TSRF.stat.full(Y,D,Cov.aug[,-(1:(Q-1-vio.index))],weight.full)
  Coef.naive[i,2] <- full.outputs$betaHat
  sd.naive[i,2] <- full.outputs$sd
  rm(weight.full)
  rm(forest.full)
  
  rm(forest.2)
  
}


Coef.matrix.basis <- na.omit(Coef.matrix.basis)
sd.matrix.basis <- na.omit(sd.matrix.basis)
apply(Coef.matrix.basis,2,mean)
apply(sd.matrix.basis,2,mean)
### coverage of fixed space for basis
Cov.mat.basis <- matrix(NA,nrow(sd.matrix.basis),ncol(sd.matrix.basis))
colnames(Cov.mat.basis) <- colnames(sd.matrix.basis)
for (j in 1:ncol(sd.matrix.basis)) {
  Cov.mat.basis[,j] <- ifelse(Coef.matrix.basis[,j]-1.96*sd.matrix.basis[,j]<=beta &
                                beta<=Coef.matrix.basis[,j]+1.96*sd.matrix.basis[,j],1,0)
}
apply(Cov.mat.basis,2,mean)



### coverage of oracle violation space
Cov.oracle = matrix(NA,nsim,ncol(sd.oracle))
colnames(Cov.oracle) = colnames(sd.oracle)
for (j in 1:ncol(sd.oracle)) {
  Cov.oracle[,j] = ifelse(Coef.oracle[,j]-1.96*sd.oracle[,j]<=beta &
                            beta<=Coef.oracle[,j]+1.96*sd.oracle[,j],1,0)
}
apply(Cov.oracle,2,mean)


### coverage of selection method
apply(Coef.robust,2,mean)
apply(sd.robust,2,mean)

Cov.robust = matrix(NA,nsim,ncol(sd.robust))
colnames(Cov.robust) = colnames(sd.robust)
for (j in 1:ncol(sd.robust)) {
  Cov.robust[,j] = ifelse(Coef.robust[,j]-1.96*sd.robust[,j]<=beta &
                            beta<=Coef.robust[,j]+1.96*sd.robust[,j],1,0)
}
apply(Cov.robust,2,mean)



### coverage of naive methods
Cov.naive <- matrix(NA,nsim,ncol(sd.naive))
colnames(Cov.naive) <- colnames(sd.naive)
for (j in 1:ncol(sd.naive)) {
  Cov.naive[,j] <- ifelse(Coef.naive[,j]-1.96*sd.naive[,j]<=beta &
                            beta<=Coef.naive[,j]+1.96*sd.naive[,j],1,0)
}
apply(Cov.naive,2,mean)


### coverage of TSLS
Cov.TSLS <- matrix(NA,nsim,ncol(sd.TSLS))
colnames(Cov.TSLS) <- colnames(sd.TSLS)
for (j in 1:ncol(sd.TSLS)) {
  Cov.TSLS[,j] <- ifelse(Coef.TSLS[,j]-1.96*sd.TSLS[,j]<=beta &
                           beta<=Coef.TSLS[,j]+1.96*sd.TSLS[,j],1,0)
}
apply(Cov.TSLS,2,mean)


###### DO NOT MODIFY ABOVE #####




