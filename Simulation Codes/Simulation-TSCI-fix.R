# simulation of fix/ee estimator


library(MASS)
source('~/Source-RF-fix.R', encoding = 'UTF-8')
source('~/Source-RF-hetero.R', encoding = 'UTF-8')


# function to generate covariance matrix
A1gen<-function(rho,p){
  A1=matrix(0,p,p)
  for(i in 1:p){
    for(j in 1:p){
      A1[i,j]<-rho^(abs(i-j))
    }
  }
  A1
}



###### dimension change across 10,20
p = 20
####please set n = 5000
n = 1000
### setting, change across 1,2
f.index = 1
##### change the interaction across 0, 0.5, 1
inter.val = 0
#### violation index, change across 1
vio.index = 1
### 
error.setting = 1
### change round across 1,2,3,4,5,6,7,8,9,10
round = 1
#### a denotes the IV strength, set as 1
a = 0.125
##### tau denotes the violation strength, set as 1
tau = 1
#### the number of simulations, make sure there is 500 in total
#### if you specify round to be 1,2,...,10, then set nsim as 50
nsim = 20


######### DO NOT MODIFY BELOW #########
f_1 <- function(x){x+a*(x^2+0.125*x^4) -25/12}
# f_2 <- function(x){exp(2*a+x+0.5*x^3)+a*(x^2+0.5*x^4)}
f_2 <- function(x){a*(1*sin(2*pi*x) + 1.5*cos(2*pi*x))}
rho1=0.5
Cov<-(A1gen(rho1,p+1))
mu<-rep(0,p+1)
beta=1
alpha<-as.matrix(rep(-0.3,p))
gamma<-as.matrix(rep(0.2,p))
inter<-as.matrix(c(rep(inter.val,5),rep(0,p-5)))


### number of violation space to consider
if (f.index==3) {
  Q = 2
} else {
  Q = 4
}


estimator.names = c(paste("RF-q",0:(Q-1),sep=""),paste("Cor-q",0:(Q-1),sep=""))
### the fixed violation space for RF
Coef.matrix.rf=matrix(NA,nrow=nsim,ncol=2*Q)
sd.matrix.rf=matrix(NA,nrow=nsim,ncol=2*Q)
colnames(Coef.matrix.rf) = colnames(sd.matrix.rf) = estimator.names


# IV strength test
iv.str = iv.thol = matrix(NA,nrow=nsim,ncol=Q)
colnames(iv.str) = colnames(iv.thol) = paste("q",0:(Q-1),sep="")


# variance of epsilon
SigmaSqY = matrix(NA, nsim, Q)
colnames(SigmaSqY) = estimator.names[1:Q]
# variance of delta
SigmaSqD = matrix(NA, nsim, 1)


### selection
Qmax  = matrix(NA,nsim,1)
colnames(Qmax) = "RF"
q.comp = q.robust = SigmaSqY.Qmax = Qmax
Coef.robust = matrix(NA,nsim,4)
sd.robust = matrix(NA,nsim,4)
colnames(Coef.robust) = colnames(sd.robust) = c("RF-comp","RF-Cor-comp","RF-robust","RF-Cor-robust")

### oracle estimator
Coef.oracle = sd.oracle = matrix(NA,nsim,2)
colnames(Coef.oracle) = colnames(sd.oracle) = c("RF","RF-Cor")

### weak IV problem
run.OLS = weak.iv = matrix(NA,nsim,1)
colnames(run.OLS) = colnames(weak.iv) = "RF"
trace.T = explained.iv = SigmaSqY



names.fix <- paste("Fix-q",0:(Q-1),sep="")
### the fixed violation space for RF
Coef.matrix.fix<-matrix(NA,nrow=nsim,ncol=Q)
sd.matrix.fix<-matrix(NA,nrow=nsim,ncol=Q)
colnames(Coef.matrix.fix) <- colnames(sd.matrix.fix) <- names.fix



# variance of epsilon
SigmaSqY.fix <- SigmaYD.fix <- D.RSS.fix = Cov.YD = matrix(NA, nsim, Q)
colnames(SigmaSqY.fix) <- colnames(SigmaYD.fix) <- colnames(D.RSS.fix) <- colnames(Cov.YD) <- names.fix[1:Q]
# variance of delta
SigmaSqD.fix = matrix(NA, nsim, 1)



for(i in 1:nsim) {
  print(i)
  #### generate the data
  mu.error<-rep(0,2)
  Cov.error<-matrix(c(1,0.5,0.5,1),2,2)
  Error<-mvrnorm(n, mu.error, Cov.error)
  W.original<-mvrnorm(n, mu, Cov)
  W<-pnorm(W.original)
  Z<-W[,1]
  if (f.index==3) {
    Z = ifelse(Z>0.6,1,0)
  } else {
    Z = 4*(Z-0.5)
  }
  X=W[,-1]
  
  ###### generate the data for the treatment variable D
  if(f.index==1){
    D=f_1(Z)+X%*%alpha+Z*X%*%inter+Error[,1]
  }
  if(f.index==2){
    D=f_2(Z)+X%*%alpha+Z*X%*%inter+Error[,1]
  }
  if(f.index==3){
    D=Z*a+X%*%alpha+Z*X%*%inter+Error[,1]
  }
  ####### generate the outcome variable
  # if(vio.index==0){
  #   Y=D*beta+ X%*% gamma+Error[,2]
  # }
  if(vio.index==1){
    Y=D*beta+ tau*Z+ X%*% gamma+Error[,2]
  }
  if(vio.index==2){
    Y=D*beta+ tau*(Z^2+Z-1)+ X%*% gamma+Error[,2] # difficult if no Z
  }
  # if(vio.index==3){
  #   Y=D*beta+ tau*Z^3+ X%*% gamma+Error[,2]
  # }
  # if(vio.index==4){
  #   Y=D*beta+ tau*Z+ Z^5/40+ X%*% gamma+Error[,2]
  # }
  #if(vio.index==5){
  #  Y=D*beta+ tau*(Z^2-1)+ 1/20*sin(Z/2)+ X%*% gamma+Error[,2]
  #}
  
  
  
  
  ### random forest based methods
  # mtry from p/3 to 2*p/3
  mtry = seq(round((p+1)/3), round(2*(p+1)/3), by=1)
  ### set max.depth and min.node.size for tuning
  ### larger max.depth and smaller min.node.size means more complex trees
  max.depth <- 0; min.node.size <- c(5,10,20)
  ### Data splitting random forest
  # use 2 to denote 2split
  forest.fix <- TSCI.RF.fit(D,Z,X,num.trees = 200,mtry=mtry,max.depth=max.depth,min.node.size=min.node.size,split.prop=2/3)
  A1.ind <- forest.fix$A1.ind
  # weight matrix
  weight.2 <- TSCI.RF.weight(forest.fix$nodes.A1)
  rm(forest.fix)
  
  Cov.aug <- X
  for (q in 1:(Q-1)) {
    Cov.aug<-cbind(Z^q,Cov.aug)
  }
  
  
  ### fix estimator
  outputs.fix <- TSCI.RF.fix(Y, D, Cov.aug, A1.ind, weight.2, Q=Q)
  ### outputs
  Coef.matrix.fix[i,] <- outputs.fix$Coef.vec
  sd.matrix.fix[i,] <- outputs.fix$sd.vec
  SigmaSqY.fix[i,] <- outputs.fix$SigmaSqY
  SigmaSqD.fix[i] <- outputs.fix$SigmaSqD
  SigmaYD.fix[i,] <- outputs.fix$SigmaYD
  D.RSS.fix[i,] <- outputs.fix$denom
  Cov.YD[i,] <- outputs.fix$numer
  
  rm(weight.2)
  
  
  
  ### TSCI
  if (f.index==3) {
    outputs.2 = TSCI.RF(Y,D,Z,X,vio.space=Z)
  } else {
    outputs.2 = TSCI.RF(Y,D,Z,X)
  }
  
  ### outputs
  Coef.matrix.rf[i,] = outputs.2$Coef.vec
  sd.matrix.rf[i,] = outputs.2$sd.vec
  SigmaSqY[i,] = outputs.2$SigmaSqY
  SigmaSqD[i] = outputs.2$SigmaSqD
  # SigmaYD.Qmax[i,] = outputs.2$SigmaYD.Qmax
  iv.str[i,] = outputs.2$iv.str; iv.thol[i,] = outputs.2$iv.thol;
  
  Coef.robust[i,] = outputs.2$Coef.robust; sd.robust[i,] = outputs.2$sd.robust;
  SigmaSqY.Qmax[i,1] = outputs.2$SigmaSqY.Qmax
  Qmax[i,1] = outputs.2$Qmax; q.comp[i,1] = outputs.2$q.comp; q.robust[i,1] = outputs.2$q.robust;
  Coef.oracle[i,1] = Coef.matrix.rf[i,vio.index+1]
  Coef.oracle[i,2] = Coef.matrix.rf[i,vio.index+Q+1]
  sd.oracle[i,1] = sd.matrix.rf[i,vio.index+1]
  sd.oracle[i,2] = sd.matrix.rf[i,vio.index+Q+1]
  run.OLS[i,1] = outputs.2$run.OLS; weak.iv[i,1] = outputs.2$weak.iv
  trace.T[i,] = outputs.2$trace.T
  explained.iv[i,] = outputs.2$explained.iv
  
  

}



apply(Coef.matrix.fix,2,mean)
apply(sd.matrix.fix,2,mean)
### coverage of fixed space for ramdom forest
Cov.mat.fix <- matrix(NA,nsim,ncol(sd.matrix.fix))
colnames(Cov.mat.fix) <- colnames(sd.matrix.fix)
for (j in 1:ncol(sd.matrix.fix)) {
  Cov.mat.fix[,j] <- ifelse(Coef.matrix.fix[,j]-1.96*sd.matrix.fix[,j]<=beta &
                             beta<=Coef.matrix.fix[,j]+1.96*sd.matrix.fix[,j],1,0)
}
apply(Cov.mat.fix,2,mean)




apply(Coef.matrix.rf,2,mean)
apply(sd.matrix.rf,2,mean)
### coverage of fixed space for ramdom forest
Cov.mat.rf <- matrix(NA,nsim,ncol(sd.matrix.rf))
colnames(Cov.mat.rf) <- colnames(sd.matrix.rf)
for (j in 1:ncol(sd.matrix.rf)) {
  Cov.mat.rf[,j] <- ifelse(Coef.matrix.rf[,j]-1.96*sd.matrix.rf[,j]<=beta &
                             beta<=Coef.matrix.rf[,j]+1.96*sd.matrix.rf[,j],1,0)
}
apply(Cov.mat.rf,2,mean)



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



######### DO NOT MODIFY ABOVE ####

# set this to your own project directory
# setwd("~/TSCI-fix/resource")
# filename <- paste("TSCç\I-fix-Setting",f.index,"-Error",error.setting,"-Interaction",inter.val,"-Violation",vio.index,"-p",p,"-n",n,"-round",round,".RData",sep="")
# save.image(filename)
# 
