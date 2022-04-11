source("Source-RF-hetero-MRule.R",encoding="UTF-8")

### Maimonides' Rule
dat5 <- read.csv("~/Desktop/projects/Real-Data-TSCI/Maimonides' Rule/cleaned data/Fifth_school_level.csv")
dat5 <- na.omit(dat5)

Y = dat5$Avgmath
X = dat5$PercentDisadvantaged
Z = dat5$classize_IV
D = dat5$Classize

# Setting 1
vio.space1 = Z
TSRF1 = TSCI.RF(Y,D,Z,X,vio.space = vio.space1,mtry = 1:2,min.node.size = seq(5,30,5))

# Setting 2
vio.space2 = dat5[,5]
TSRF2 = TSCI.RF(Y,D,Z,X,vio.space = vio.space2,mtry = 1:2,min.node.size = seq(5,30,5))


# Setting 3
vio.space2 = dat5[,c(9,5)]
TSRF3 = TSCI.RF(Y,D,Z,X,vio.space = vio.space2,mtry = 1:3,min.node.size = seq(5,30,5))


# Setting 4
vio.space2 = dat5[,10]
TSRF4 = TSCI.RF(Y,D,Z,X,vio.space = vio.space2,mtry = 1:2,min.node.size = seq(5,30,5))



attach(dat5)
library(ivmodel)
iv1 <- ivmodel(Y=Avgmath,D=Classize,Z=classize_IV,X=PercentDisadvantaged)
summary(iv1)

iv2 <- ivmodel(Y=Avgmath,D=Classize,Z=classize_IV,X=cbind(PercentDisadvantaged,Enrollment))
summary(iv2)

iv3 <- ivmodel(Y=Avgmath,D=Classize,Z=classize_IV,X=cbind(PercentDisadvantaged,Enrollment,Enrollment2))
summary(iv3)

iv4 <- ivmodel(Y=Avgmath,D=Classize,Z=classize_IV,X=cbind(PiecewiseLinearTrend))
summary(iv4)