library(kableExtra)
library(tidyverse)
# 
out.table <- data.frame(matrix(0,18,14))
colnames(out.table) <- c("vio","inter","n",
                         "Orac","Comp","Robust", # Corrected RF
                         "Orac","Comp","Robust", #Basis
                         "Regular", # TSLS
                         "Split", "Plug", "Full", # vio space
                         "TSCI-RF" # RF
)

out.table[,1] = rep(c(1,2),each = 9)
out.table[,2] = rep(c(0,0,0,0.5,0.5,0.5,1,1,1),2) # vio.index
out.table[,3] = rep(c(1000,3000,5000),6)


f.index = 1
p = 20
r = 1
for (vio.index in 1:2) {
  for (inter.val in c(0,0.5,1)) {
    for (n in c(1000,3000,5000)) {
      filename <- paste("TSCI-homo-Setting",f.index,"-Error1","-Interaction",inter.val,"-Violation",vio.index,"-p",p,"-n",n,".RData",sep="")
      load(paste("~/Desktop/TSCI-test/homo-merged/",filename,sep=""))
      out.table[r,4] <- 3.92*apply(sd.oracle.all,2,mean)[2]
      out.table[r,5] <- 3.92*apply(sd.robust.all,2,mean)[2]
      out.table[r,6] <- 3.92*apply(sd.robust.all,2,mean)[4]
      out.table[r,7] <- 3.92*apply(sd.basis.oracle.all,2,mean)
      out.table[r,8:9] <- 3.92*apply(sd.basis.robust.all,2,mean)
      out.table[r,10] = 3.92*apply(sd.TSLS.all,2,mean)
      out.table[r,11] = 3.92*apply(sd.oracle.all,2,mean)[1]
      out.table[r,12:13] = 3.92*apply(sd.naive.all,2,mean)
      out.table[r,14] = mean(q.comp.all>=1)
      r <- r + 1
    }
  }
}


# 456
kbl(out.table,"latex",caption = "", align = "c", digits = 2) %>%
  collapse_rows(columns = 1:3, latex_hline = "major") %>%
  add_header_above(c(rep(" ", 3), "TSCI-RF"=3, "Basis"=3, "TSLS", "Other RF"=3, "Invalidity")) %>%
  kable_styling(latex_options = c("scale_down"))



# fix

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
#     out.table[r,3] <- 3.92*apply(sd.oracle.all,2,mean)[2]
#     out.table[r,4] <- 3.92*apply(sd.robust.all,2,mean)[2]
#     out.table[r,5] <- 3.92*apply(sd.robust.all,2,mean)[4]
#     
#     
#     out.table[r,6] <- 3.92*apply(sd.oracle.all,2,mean)[1]
#     
#     
#     out.table[r,7:10] <- 3.92*apply(sd.matrix.fix.all,2,mean)
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