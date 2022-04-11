library(kableExtra)
library(tidyverse)

out.table <- data.frame(matrix(0,18,15))
colnames(out.table) <- c("vio","inter","n",
                         "Orac","Comp","Robust", # Corrected RF, bias
                         "Orac","Comp","Robust", # Corrected RF, Len
                         "Orac","Comp","Robust", # Corrected RF, Cover
                         "", # Invalidity
                         "Orac","Orac" # RF-Init
)

out.table[,1] = rep(c(1,2),each = 9)
out.table[,2] = rep(c(0,0,0,0.5,0.5,0.5,1,1,1),2) # vio.index
out.table[,3] = rep(c(1000,3000,5000),6)


f.index = 2
p = 20
r = 1
for (vio.index in 1:2) {
  for (inter.val in c(0,0.5,1)) {
    for (n in c(1000,3000,5000)) {
      filename <- paste("TSCI-hetero-Setting",f.index,"-Error3","-Interaction",inter.val,"-Violation",vio.index,"-p",p,"-n",n,".RData",sep="")
      load(paste("~/Desktop/TSCI-test/Error3/",filename,sep=""))
      out.table[r,4] <- abs(apply(Coef.oracle.all,2,mean)[2]-1)
      out.table[r,5] <- abs(apply(Coef.robust.all,2,mean)[2]-1)
      out.table[r,6] <- abs(apply(Coef.robust.all,2,mean)[4]-1)
      
      out.table[r,7] <- 3.92*apply(sd.oracle.all,2,mean)[2]
      out.table[r,8] <- 3.92*apply(sd.robust.all,2,mean)[2]
      out.table[r,9] <- 3.92*apply(sd.robust.all,2,mean)[4]
      
      out.table[r,10] <- apply(Cov.oracle.all,2,mean)[2]
      out.table[r,11] <- apply(Cov.robust.all,2,mean)[2]
      out.table[r,12] <- apply(Cov.robust.all,2,mean)[4]
      
      out.table[r,13] = mean(q.comp.all>=1)
      
      out.table[r,14] <- abs(apply(Coef.oracle.all,2,mean)[1]-1)
      out.table[r,15] <- apply(Cov.oracle.all,2,mean)[1]
    
      r <- r + 1
    }
  }
}


# 456
kbl(out.table,"latex",caption = "", align = "c", digits = 2) %>%
  collapse_rows(columns = 1:3, latex_hline = "major") %>%
  add_header_above(c(rep(" ", 3), "Bias"=3, "Length"=3, "Coverage"=3, "Invalidity", "Bias", "Coverage")) %>%
  add_header_above(c(rep("",3), "TSCI-RF"= 10, "RF-Init"=2)) %>%
  kable_styling(latex_options = c("scale_down"))