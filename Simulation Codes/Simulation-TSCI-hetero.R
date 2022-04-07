# simulation of TSCI in heteroskedastic settings


library(MASS)
source('~/Source-RF-hetero.R', encoding = 'UTF-8')


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
####please set n = 3000
n = 3000
### setting, change across 1,2,3
f.index = 3
##### change the interaction across 0, 0.5, 1
inter.val = 0.5
#### violation index, change across 1,2
vio.index = 1
# error distribution, 2 or 3
error.setting = 3
### change round across 1,2,3,4,5,6,7,8,9,10
round = 1
#### a denotes the IV strength, set as 1
a = 1
##### tau denotes the violation strength, set as 1
tau = 1
#### the number of simulations, make sure there is 500 in total
#### if you specify round to be 1,2,...,10, then set nsim as 50
nsim = 2


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
if (f.index==3) {
  Q = 2
} else {
  Q = 4
}
# estimator names
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
# covariance of epsilon and delta
SigmaSqY.Qmax = matrix(NA, nsim, 1)


### selection
Qmax  = matrix(NA,nsim,1)
colnames(Qmax) = "RF"
q.comp = q.robust = Qmax
Coef.robust = matrix(NA,nsim,4)
sd.robust = matrix(NA,nsim,4)
colnames(Coef.robust) = colnames(sd.robust) = c("RF-comp","RF-Cor-comp","RF-robust","RF-Cor-robust")

### oracle estimator
Coef.oracle = sd.oracle = matrix(NA,nsim,2)
colnames(Coef.oracle) = colnames(sd.oracle) = c("RF","RF-Cor")

### weak IV problem
run.OLS = weak.iv = matrix(NA,nsim,1)
colnames(run.OLS) = colnames(weak.iv) = "RF"
trace.T = explained.iv = iv.str


for(i in 1:nsim) {
  print(i)
  #### generate the data
  mu.error=rep(0,2)
  Cov.error=matrix(c(1,0.5,0.5,1),2,2)
  Error = matrix(NA,n,2)
  
  
  W.original=mvrnorm(n, mu, Cov)
  W=pnorm(W.original)
  Z=W[,1]
  if (f.index==3) {
    Z = ifelse(Z>0.6,1,0)
  } else {
    Z = 4*(Z-0.5)
  }
  X=W[,-1]
  
  
  # heteroscedastic errors
  w1 = rep(0,n)
  for (j in 1:n) {
    w1[j] = rnorm(1,0,sqrt(Z[j]^2+0.25))
  }
  w2 = rnorm(n)
  if (error.setting == 2) {
    Error[,1] = rnorm(n)
    for (j in 1:n) {
      Error[j,2] = 0.3*Error[j,1] + sqrt((1-0.3^2)/(0.86^4+1.38072^2))*(w1[j]*1.38072+w2[j]*0.86^2) 
    }
  }
  if (error.setting == 3) {
    for (j in 1:n) {
      Error[j,1] = rnorm(1,0,sqrt(Z[j]^2+0.25))
      Error[j,2] = 0.6*Error[j,1] + sqrt((1-0.6^2)/(0.86^4+1.38072^2))*(w1[j]*1.38072+w2[j]*0.86^2) 
    }
  }

  
  
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
    Y=D*beta+ tau*Z+ X%*%gamma+Error[,2]
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
  if (f.index==3) {
    outputs.2 = TSCI.RF(Y,D,Z,X,vio.space=Z)
  } else {
    outputs.2 = TSCI.RF(Y,D,Z,X)
  }
  
  
  ### outputs
  Coef.matrix.rf[i,] = outputs.2$Coef.vec
  sd.matrix.rf[i,] = outputs.2$sd.vec
  iv.str[i,] = outputs.2$iv.str; iv.thol[i,] = outputs.2$iv.thol;
  
  Coef.robust[i,] = outputs.2$Coef.robust; sd.robust[i,] = outputs.2$sd.robust;
  Qmax[i,1] = outputs.2$Qmax; q.comp[i,1] = outputs.2$q.comp; q.robust[i,1] = outputs.2$q.robust;
  Coef.oracle[i,1] = Coef.matrix.rf[i,vio.index+1]
  Coef.oracle[i,2] = Coef.matrix.rf[i,vio.index+Q+1]
  sd.oracle[i,1] = sd.matrix.rf[i,vio.index+1]
  sd.oracle[i,2] = sd.matrix.rf[i,vio.index+Q+1]
  run.OLS[i,1] = outputs.2$run.OLS; weak.iv[i,1] = outputs.2$weak.iv
  trace.T[i,] = outputs.2$trace.T
  explained.iv[i,] = outputs.2$explained.iv
  
  SigmaSqY[i,] = outputs.2$SigmaSqY
  SigmaSqD[i] = outputs.2$SigmaSqD
  SigmaSqY.Qmax[i,1] = outputs.2$SigmaSqY.Qmax
}


apply(Coef.matrix.rf,2,mean)
apply(sd.matrix.rf,2,mean)
### coverage of fixed space for ramdom forest
Cov.mat.rf = matrix(NA,nsim,ncol(sd.matrix.rf))
colnames(Cov.mat.rf) = colnames(sd.matrix.rf)
for (j in 1:ncol(sd.matrix.rf)) {
  Cov.mat.rf[,j] = ifelse(Coef.matrix.rf[,j]-1.96*sd.matrix.rf[,j]<=beta &
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
filename = paste("TSCI-hetero-Setting",f.index,"-Error",error.setting,"-Interaction",inter.val,"-Violation",vio.index,"-p",p,"-n",n,"-round",round,".RData",sep="")
save.image(filename)



