library(kableExtra)
library(tidyverse)

out.table <- data.frame(matrix(0,18,11))
colnames(out.table) <- c("Inter","vio","n",
                         "Orac","Comp","Robust", # Corrected RF
                         "Orac", # RF
                         "q0","q1","q2","q3" # vio space
)


out.table[,1] = c(rep(0,6),rep(0.5,6),rep(1,6))
out.table[,2] = rep(c(1,1,1,2,2,2),3) # vio.index
out.table[,3] = rep(c(1000,3000,5000),6)


error.setting = 3
f.index = 2
p = 20
r <- 1

for (inter.val in c(0,0.5,1)) {
  for (vio.index in 1:2) {
    for (n in c(1000,3000,5000)) {
      filename <- paste("TSCI-hetero-Setting",f.index,"-Error",error.setting,"-Interaction",inter.val,"-Violation",vio.index,"-p",p,"-n",n,".RData",sep="")
      load(paste("~/Desktop/TSCI-test/Error3/",filename,sep=""))
      out.table[r,4] <- abs(apply(Coef.oracle.all,2,mean)[2]-1)
      out.table[r,5] <- abs(apply(Coef.robust.all,2,mean)[2]-1)
      out.table[r,6] <- abs(apply(Coef.robust.all,2,mean)[4]-1)
      out.table[r,7] <- abs(apply(Coef.oracle.all,2,mean)[1]-1)
      out.table[r,8:11] <- abs(apply(Coef.matrix.rf.all,2,mean)[5:8]-1)
      r <- r + 1
    }
  }
}


# 456
kbl(out.table,"latex",caption = "", align = "c", digits = 2) %>%
  collapse_rows(columns = 1:3, latex_hline = "major") %>%
  add_header_above(c(rep(" ", 3), "RF-Cor"=3, "RF", "Violation Space"=4)) %>%
  kable_styling(latex_options = c("scale_down"))
