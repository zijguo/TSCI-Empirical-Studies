library(kableExtra)
library(tidyverse)





out.table <- data.frame(matrix(0,15,12))
colnames(out.table) <- c("a","n",
                         "TSCI-RF","RF-Init","RF-Fix", # bias
                         "TSCI-RF","RF-Init","RF-Fix", # length
                         "TSCI-RF","RF-Fix", # coverage
                         "TSCI-RF", "TSCI-RF"
)

out.table[,1] = rep(c(0.1,0.12,0.125,0.15,0.2),each = 3)
out.table[,2] = rep(c(1000,3000,5000),5)



f.index = 1
inter.val = 0
vio.index = 1
p = 20
r = 1

for (a in c(0.1,0.12,0.125,0.15,0.2)) {
  for (n in c(1000,3000,5000)) {
    filename <- paste("TSCI-fix-Setting",f.index,"-Str",a,"-Error1","-Interaction",inter.val,"-Violation",vio.index,"-p",p,"-n",n,".RData",sep="")
    load(paste("~/Desktop/TSCI-test/fix/",filename,sep=""))
    out.table[r,3] <- abs(apply(Coef.oracle.all,2,mean)[2]-1)
    out.table[r,4] <- abs(apply(Coef.oracle.all,2,mean)[1]-1)
    out.table[r,5] <- abs(apply(Coef.matrix.fix.all,2,mean)[2]-1)



    out.table[r,6] <- apply(Cov.oracle.all,2,mean)[2]
    out.table[r,7] <- apply(Cov.oracle.all,2,mean)[1]
    out.table[r,8] <- apply(Cov.mat.fix.all,2,mean)[2]

    out.table[r,9] <- 3.92*apply(sd.oracle.all,2,mean)[2]
    out.table[r,10] <- 3.92*apply(sd.matrix.fix.all,2,mean)[2]


    out.table[r,11] = apply(iv.str.all,2,mean)[2]
    out.table[r,12] = mean(run.OLS.all) + mean(weak.iv.all)



    r <- r + 1
  }
}


# 456
kbl(out.table,"latex",caption = "", align = "c", digits = 2) %>%
  collapse_rows(columns = 1:2, latex_hline = "major") %>%
  add_header_above(c(rep(" ", 2), "Bias"=3, "Coverage"=3, "CI Length"=2, "IV Strength", "Weak IV")) %>%
  kable_styling(latex_options = c("scale_down"))




# library(kableExtra)
# library(tidyverse)
# 
# 
# # fix and RF
# out.table <- data.frame(matrix(0,18,13))
# colnames(out.table) <- c("vio","inter","n",
#                          "TSCI-RF","RF-Init","RF-Fix", # bias
#                          "TSCI-RF","RF-Init","RF-Fix", # coverage
#                          "TSCI-RF","RF-Fix", # length
#                          "TSCI-RF", "TSCI-RF"
# )
# 
# out.table[,1] = rep(c(1,2),each = 9)
# out.table[,2] = rep(c(0,0,0,0.5,0.5,0.5,1,1,1),2) # vio.index
# out.table[,3] = rep(c(1000,3000,5000),6)
# 
# 
# f.index = 1
# p = 20
# r = 1
# for (vio.index in 1:2) {
#   for (inter.val in c(0,0.5,1)) {
#     for (n in c(1000,3000,5000)) {
#       filename <- paste("TSCI-homo-Setting",f.index,"-Error1","-Interaction",inter.val,"-Violation",vio.index,"-p",p,"-n",n,".RData",sep="")
#       load(paste("~/Desktop/TSCI-test/fix2&RF-Setting1/",filename,sep=""))
#       out.table[r,4] <- abs(apply(Coef.oracle.all,2,mean)[2]-beta)
#       out.table[r,5] <- abs(apply(Coef.oracle.all,2,mean)[1]-beta)
#       out.table[r,6] <- abs(apply(Coef.matrix.fix.all,2,mean)[vio.index+1]-beta)
# 
# 
# 
#       out.table[r,7] <- apply(Cov.oracle.all,2,mean)[2]
#       out.table[r,8] <- apply(Cov.oracle.all,2,mean)[1]
#       out.table[r,9] <- apply(Cov.mat.fix.all,2,mean)[vio.index+1]
# 
#       out.table[r,10] <- 3.92*apply(sd.oracle.all,2,mean)[2]
#       out.table[r,11] <- 3.92*apply(sd.matrix.fix.all,2,mean)[vio.index+1]
# 
# 
#       out.table[r,12] = apply(iv.str.all,2,mean)[vio.index+1]
#       out.table[r,13] = mean(run.OLS.all) + mean(weak.iv.all)
# 
# 
#       r <- r + 1
#     }
#   }
# }
# 
# 
# # 456
# kbl(out.table,"latex",caption = "", align = "c", digits = 2) %>%
#   collapse_rows(columns = 1:3, latex_hline = "major") %>%
#   add_header_above(c(rep(" ", 3), "Bias"=3, "Coverage"=3, "CI Length"=2, "IV Strength", "Weak IV")) %>%
#   kable_styling(latex_options = c("scale_down"))





# # fix and RF binary
# 
# 
# library(kableExtra)
# library(tidyverse)
# 
# 
# # fix and RF
# out.table <- data.frame(matrix(0,9,12))
# colnames(out.table) <- c("inter","n",
#                          "TSCI-RF","RF-Init","RF-Fix", # bias
#                          "TSCI-RF","RF-Init","RF-Fix", # coverage
#                          "TSCI-RF","RF-Fix", # length
#                          "TSCI-RF", "TSCI-RF"
# )
# 
# 
# out.table[,2] = c(0,0,0,0.5,0.5,0.5,1,1,1) # vio.index
# out.table[,3] = rep(c(1000,3000,5000),3)
# 
# 
# f.index = 3
# p = 20
# r = 1
# for (vio.index in 1) {
#   for (inter.val in c(0,0.5,1)) {
#     for (n in c(1000,3000,5000)) {
#       filename <- paste("TSCI-homo-Setting",f.index,"-Error1","-Interaction",inter.val,"-Violation",vio.index,"-p",p,"-n",n,".RData",sep="")
#       load(paste("~/Desktop/TSCI-test/fix2&RF-binary/",filename,sep=""))
#       out.table[r,3] <- abs(apply(Coef.oracle.all,2,mean)[2]-beta)
#       out.table[r,4] <- abs(apply(Coef.oracle.all,2,mean)[1]-beta)
#       out.table[r,5] <- abs(apply(Coef.matrix.fix.all,2,mean)[vio.index+1]-beta)
#       
#       
#       
#       out.table[r,6] <- apply(Cov.oracle.all,2,mean)[2]
#       out.table[r,7] <- apply(Cov.oracle.all,2,mean)[1]
#       out.table[r,8] <- apply(Cov.mat.fix.all,2,mean)[vio.index+1]
#       
#       out.table[r,9] <- 3.92*apply(sd.oracle.all,2,mean)[2]
#       out.table[r,10] <- 3.92*apply(sd.matrix.fix.all,2,mean)[vio.index+1]
#       
#       
#       out.table[r,12] = apply(iv.str.all,2,mean)[vio.index+1]
#       out.table[r,12] = mean(run.OLS.all) + mean(weak.iv.all)
#       
#       
#       r <- r + 1
#     }
#   }
# }
# 
# 
# # 456
# kbl(out.table,"latex",caption = "", align = "c", digits = 3) %>%
#   collapse_rows(columns = 1:2, latex_hline = "major") %>%
#   add_header_above(c(rep(" ", 2), "Bias"=3, "Coverage"=3, "CI Length"=2, "IV Strength", "Weak IV")) %>%
#   kable_styling(latex_options = c("scale_down"))









