# library(kableExtra)
# library(tidyverse)
# 
# out.table <- data.frame(matrix(0,18,14))
# colnames(out.table) <- c("vio","inter","n",
#                          "Orac","Comp","Robust", # Corrected RF
#                          "Orac","Comp","Robust", #Basis
#                          "Regular", # TSLS
#                          "Split", "Plug", "Full", # vio space
#                          "TSCI-RF" # RF
# )
# 
# out.table[,1] = rep(c(1,2),each = 9)
# out.table[,2] = rep(c(0,0,0,0.5,0.5,0.5,1,1,1),2) # vio.index
# out.table[,3] = rep(c(1000,3000,5000),6)
# 
# 
# f.index = 2
# p = 20
# r = 1
# for (vio.index in 1:2) {
#   for (inter.val in c(0,0.5,1)) {
#     for (n in c(1000,3000,5000)) {
#       filename <- paste("TSCI-homo-Setting",f.index,"-Error1","-Interaction",inter.val,"-Violation",vio.index,"-p",p,"-n",n,".RData",sep="")
#       load(paste("~/Desktop/TSCI-test/homo-merged/",filename,sep=""))
#       out.table[r,4] <- abs(apply(Coef.oracle.all,2,mean)[2]-1)
#       out.table[r,5] <- abs(apply(Coef.robust.all,2,mean)[2]-1)
#       out.table[r,6] <- abs(apply(Coef.robust.all,2,mean)[4]-1)
#       out.table[r,7] <- abs(apply(Coef.basis.oracle.all,2,mean)-1)
#       out.table[r,8:9] <- abs(apply(Coef.basis.robust.all,2,mean)-1)
#       out.table[r,10] = abs(apply(Coef.TSLS.all,2,mean)-1)
#       out.table[r,11] = abs(apply(Coef.oracle.all,2,mean)[1]-1)
#       out.table[r,12:13] = abs(apply(Coef.naive.all,2,mean)-1)
#       out.table[r,14] = mean(q.comp.all>=1)
#       r <- r + 1
#     }
#   }
# }
# 
# 
# # 456
# kbl(out.table,"latex",caption = "", align = "c", digits = 2) %>%
#   collapse_rows(columns = 1:3, latex_hline = "major") %>%
#   add_header_above(c(rep(" ", 3), "TSCI-RF"=3, "Basis"=3, "TSLS", "Other RF"=3, "Invalidity")) %>%
#   kable_styling(latex_options = c("scale_down"))





# ### Binary Treatment
# 
# library(kableExtra)
# library(tidyverse)
# 
# out.table <- data.frame(matrix(0,18,16))
# colnames(out.table) <- c("Setting","Inter","n",
#                          "Orac","Comp","Robust", # Corrected RF
#                          "Orac","Comp","Robust", # Corrected RF
#                          "Orac","Comp","Robust", # Corrected RF
#                          "TSCI-RF",
#                          "Orac","Orac", # RF
#                          "" # weak iv
# 
# )
# 
# out.table[,1] = rep(c(1,2),each = 9)
# out.table[,2] = rep(c(0,0,0,0.5,0.5,0.5,1,1,1),2) # vio.index
# out.table[,3] = rep(c(1000,3000,5000),6)
# 
# 
# vio.index = 1
# p = 20
# r = 1
# for (f.index in 1:2) {
#   for (inter.val in c(0,0.5,1)) {
#     for (n in c(1000,3000,5000)) {
#       filename <- paste("TSCI-BiCon-Setting",f.index,"-Error1","-Interaction",inter.val,"-Violation",vio.index,"-p",p,"-n",n,".RData",sep="")
#       load(paste("~/Desktop/TSCI-test/BiCon/",filename,sep=""))
#       out.table[r,4] <- abs(apply(Coef.oracle.all,2,mean)[2]-1)
#       out.table[r,5] <- abs(apply(Coef.robust.all,2,mean)[2]-1)
#       out.table[r,6] <- abs(apply(Coef.robust.all,2,mean)[4]-1)
# 
#       out.table[r,7] <- 3.92*apply(sd.oracle.all,2,mean)[2]
#       out.table[r,8] <- 3.92*apply(sd.robust.all,2,mean)[2]
#       out.table[r,9] <- 3.92*apply(sd.robust.all,2,mean)[4]
# 
#       out.table[r,10] <- apply(Cov.oracle.all,2,mean)[2]
#       out.table[r,11] <- apply(Cov.robust.all,2,mean)[2]
#       out.table[r,12] <- apply(Cov.robust.all,2,mean)[4]
# 
#       out.table[r,13] = mean(q.comp.all>=1)
# 
#       out.table[r,14] <- abs(apply(Coef.oracle.all,2,mean)[1]-1)
#       out.table[r,15] <- apply(Cov.oracle.all,2,mean)[1]
#       out.table[r,16] <- mean(weak.iv.all) + mean(run.OLS.all)
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
#   add_header_above(c(rep(" ", 3), "Bias"=3, "Length"=3, "Coverage"=3, "Invalidity", "Bias", "Coverage", "Weak IV")) %>%
#   add_header_above(c(rep(" ", 3), "TSCI-RF"=10,"RF-Init"=2,"TSCI-RF")) %>%
#   kable_styling(latex_options = c("scale_down"))






# ### Binary Treatment, larger interaction
# 
# library(kableExtra)
# library(tidyverse)
# 
# out.table <- data.frame(matrix(0,12,16))
# colnames(out.table) <- c("Setting","Inter","n",
#                          "Orac","Comp","Robust", # Corrected RF
#                          "Orac","Comp","Robust", # Corrected RF
#                          "Orac","Comp","Robust", # Corrected RF
#                          "TSCI-RF",
#                          "Orac","Orac", # RF
#                          "" # weak iv
#                          
# )
# 
# out.table[,1] = rep(c(1,2),each = 6)
# out.table[,2] = rep(c(1.5,1.5,1.5,2,2,2),2) # vio.index
# out.table[,3] = rep(c(1000,3000,5000),4)
# 
# 
# vio.index = 1
# p = 20
# r = 1
# for (f.index in 1:2) {
#   for (inter.val in c(1.5,2)) {
#     for (n in c(1000,3000,5000)) {
#       filename <- paste("TSCI-BiCon-Setting",f.index,"-Error1","-Interaction",inter.val,"-Violation",vio.index,"-p",p,"-n",n,".RData",sep="")
#       load(paste("~/Desktop/TSCI-test/BiCon-largeInter/",filename,sep=""))
#       out.table[r,4] <- abs(apply(Coef.oracle.all,2,mean)[2]-1)
#       out.table[r,5] <- abs(apply(Coef.robust.all,2,mean)[2]-1)
#       out.table[r,6] <- abs(apply(Coef.robust.all,2,mean)[4]-1)
#       
#       out.table[r,7] <- 3.92*apply(sd.oracle.all,2,mean)[2]
#       out.table[r,8] <- 3.92*apply(sd.robust.all,2,mean)[2]
#       out.table[r,9] <- 3.92*apply(sd.robust.all,2,mean)[4]
#       
#       out.table[r,10] <- apply(Cov.oracle.all,2,mean)[2]
#       out.table[r,11] <- apply(Cov.robust.all,2,mean)[2]
#       out.table[r,12] <- apply(Cov.robust.all,2,mean)[4]
#       
#       out.table[r,13] = mean(q.comp.all>=1)
#       
#       out.table[r,14] <- abs(apply(Coef.oracle.all,2,mean)[1]-1)
#       out.table[r,15] <- apply(Cov.oracle.all,2,mean)[1]
#       out.table[r,16] <- mean(weak.iv.all) + mean(run.OLS.all)
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
#   add_header_above(c(rep(" ", 3), "Bias"=3, "Length"=3, "Coverage"=3, "Invalidity", "Bias", "Coverage", "Weak IV")) %>%
#   add_header_above(c(rep(" ", 3), "TSCI-RF"=10,"RF-Init"=2,"TSCI-RF")) %>%
#   kable_styling(latex_options = c("scale_down"))






# # fix
# out.table <- data.frame(matrix(0,15,10))
# colnames(out.table) <- c("Strength","n",
#                          "Orac","Comp","Robust", # Corrected RF
#                          "Orac",
#                          "q0","q1","q2","q3"
# )
# 
# out.table[,1] = rep(c(0.1,0.12,0.125,0.15,0.2),each = 3)
# out.table[,2] = rep(c(1000,3000,5000),5)
# 
# 
# 
# f.index = 1
# inter.val = 0
# vio.index = 1
# p = 20
# r = 1
# 
# for (a in c(0.1,0.12,0.125,0.15,0.2)) {
#   for (n in c(1000,3000,5000)) {
#     filename <- paste("TSCI-fix-Setting",f.index,"-Str",a,"-Error1","-Interaction",inter.val,"-Violation",vio.index,"-p",p,"-n",n,".RData",sep="")
#     load(paste("~/Desktop/TSCI-test/fix/",filename,sep=""))
#     out.table[r,3] <- abs(apply(Coef.oracle.all,2,mean)[2]-1)
#     out.table[r,4] <- abs(apply(Coef.robust.all,2,mean)[2]-1)
#     out.table[r,5] <- abs(apply(Coef.robust.all,2,mean)[4]-1)
# 
# 
#     out.table[r,6] <- abs(apply(Coef.oracle.all,2,mean)[1]-1)
# 
# 
#     out.table[r,7:10] <- abs(apply(Coef.matrix.fix.all,2,mean)-1)
# 
# 
# 
#     r <- r + 1
#   }
# }
# 
# 
# # 456
# kbl(out.table,"latex",caption = "", align = "c", digits = 3) %>%
#   collapse_rows(columns = 1:2, latex_hline = "major") %>%
#   add_header_above(c(rep(" ", 2), "TSCI-RF"=3, "RF-Init", "RF-Fix"=4)) %>%
#   kable_styling(latex_options = c("scale_down"))





















