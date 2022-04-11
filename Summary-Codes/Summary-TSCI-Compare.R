# ### Continuous Treatment
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
#                          "","","" # RF
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
#       filename <- paste("TSCI-Compare-Setting",f.index,"-Error4","-Interaction",inter.val,"-Violation",vio.index,"-p",p,"-n",n,".RData",sep="")
#       load(paste("~/Desktop/TSCI-test/TSCI-Compare/",filename,sep=""))
#       out.table[r,4] <- abs(apply(Coef.oracle.all,2,mean)[2]-beta)
#       out.table[r,5] <- abs(apply(Coef.robust.all,2,mean)[2]-beta)
#       out.table[r,6] <- abs(apply(Coef.robust.all,2,mean)[4]-beta)
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
#       out.table[r,14] <- abs(apply(Coef.genius.all,2,mean)-beta)
#       out.table[r,15] = 3.92*apply(sd.genius.all,2,mean)
#       out.table[r,16] <- apply(Cov.genius.all,2,mean)
# 
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
#   add_header_above(c(rep(" ", 3), "Bias"=3, "Length"=3, "Coverage"=3, "Invalidity", "Bias", "Length", "Coverage")) %>%
#   add_header_above(c(rep(" ", 3), "TSCI-RF"=10,"MR Genius"=3)) %>%
#   kable_styling(latex_options = c("scale_down"))





# Binary Treatment
### Continuous Treatment

library(kableExtra)
library(tidyverse)

out.table <- data.frame(matrix(0,9,16))
colnames(out.table) <- c("Inter","n",
                         "Orac","Comp","Robust", # Corrected RF
                         "Orac","Comp","Robust", # Corrected RF
                         "Orac","Comp","Robust", # Corrected RF
                         "TSCI-RF",
                         "","","" # RF
                         
)


out.table[,1] = c(0,0,0,0.5,0.5,0.5,1,1,1) # vio.index
out.table[,2] = rep(c(1000,3000,5000),3)


vio.index = 1
p = 20
r = 1

for (inter.val in c(0,0.5,1)) {
  for (n in c(1000,3000,5000)) {
    filename <- paste("TSCI-Compare-BiTreat","-Error1","-Interaction",inter.val,"-Violation",vio.index,"-p",p,"-n",n,".RData",sep="")
    load(paste("~/Desktop/TSCI-test/TSCI-Compare-BiTreat/",filename,sep=""))
    out.table[r,3] <- abs(apply(Coef.oracle.all,2,mean)[2]-beta)
    out.table[r,4] <- abs(apply(Coef.robust.all,2,mean)[2]-beta)
    out.table[r,5] <- abs(apply(Coef.robust.all,2,mean)[4]-beta)
    
    out.table[r,6] <- 3.92*apply(sd.oracle.all,2,mean)[2]
    out.table[r,7] <- 3.92*apply(sd.robust.all,2,mean)[2]
    out.table[r,8] <- 3.92*apply(sd.robust.all,2,mean)[4]
    
    out.table[r,9] <- apply(Cov.oracle.all,2,mean)[2]
    out.table[r,10] <- apply(Cov.robust.all,2,mean)[2]
    out.table[r,11] <- apply(Cov.robust.all,2,mean)[4]
    
    out.table[r,12] = mean(q.comp.all>=1)
    
    out.table[r,13] <- abs(apply(Coef.genius.all,2,mean)-beta)
    out.table[r,14] = 3.92*apply(sd.genius.all,2,mean)
    out.table[r,15] <- apply(Cov.genius.all,2,mean)
    
    out.table[r,16] = mean(weak.iv.all) + mean(run.OLS.all)
    
    
    r <- r + 1
  }
}



# 456
kbl(out.table,"latex",caption = "", align = "c", digits = 2) %>%
  collapse_rows(columns = 1:2, latex_hline = "major") %>%
  add_header_above(c(rep(" ", 2), "Bias"=3, "Length"=3, "Coverage"=3, "Invalidity", "Bias", "Length", "Coverage")) %>%
  add_header_above(c(rep(" ", 2), "TSCI-RF"=10,"MR Genius"=3)) %>%
  kable_styling(latex_options = c("scale_down"))











