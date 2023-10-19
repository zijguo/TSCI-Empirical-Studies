library(xtable)

# table -------------------------------------------------------------------

ana_tab = NULL
for (vio in c(1,2)) {
    for (a in c(0,0.5,1)) {
        for (n in c(1000,3000,5000)) {
            total = NULL
            for (round in 1:25) {
                filename = paste0("results-model1/results-Simulation-TSCI-invalidIV-model1"
                                  , "-a", a
                                  , "-vio", vio
                                  , "-n", n
                                  , "-round", round, ".rds")
                one_round = readRDS(filename)
                total = rbind(total, one_round)
            }
            print(dim(total))
            beta0 = 1
            bias_oracle = bias_comp = bias_robust = bias_tsls = bias_init = NULL
            bias_plug_hetero = bias_full_hetero = NULL
            qcomp = invalidity = NULL
            for (i in 1:nrow(total)) {
                bias_oracle = c(bias_oracle, total[i,"res_tsci"][[1]]$Coef.vec[5+vio] - beta0)
                if (total[i,"res_tsci"][[1]]$weak.iv == T) {
                    qcomp = c(qcomp, 0)
                    qrobust = 0
                    invalidity = c(invalidity, 0)
                } else {
                    qcomp = c(qcomp, total[i,"res_tsci"][[1]]$q.comp)
                    qrobust = total[i,"res_tsci"][[1]]$q.robust
                    invalidity = c(invalidity, total[i,"res_tsci"][[1]]$invalidity)
                }
                bias_comp = c(bias_comp, total[i,"res_tsci"][[1]]$Coef.vec[5+qcomp[i]] - beta0)
                bias_robust = c(bias_robust, total[i,"res_tsci"][[1]]$Coef.vec[5+qrobust] - beta0)
                bias_tsls = c(bias_tsls, total[i,"res_tsls"][[1]][1] - beta0)
                bias_init = c(bias_init, total[i,"res_tsci"][[1]]$Coef.vec[vio+1] - beta0)
                bias_plug_hetero = c(bias_plug_hetero, total[i,"res_plug_hetero"][[1]]$betaHat - beta0)
                bias_full_hetero = c(bias_full_hetero, total[i,"res_full_hetero"][[1]]$betaHat - beta0)
            }
            record = c(
                vio = vio
                , a = a
                , n = n
                , abs(mean(bias_oracle))
                , abs(mean(bias_comp))
                , abs(mean(bias_robust))
                , mean(invalidity)
                , sum(qcomp==0) / nrow(total)
                , sum(qcomp==1) / nrow(total)
                , sum(qcomp==2) / nrow(total)
                , sum(qcomp==3) / nrow(total)
                , abs(mean(bias_tsls))
                , abs(mean(bias_init))
                , abs(mean(bias_plug_hetero))
                , abs(mean(bias_full_hetero))
            )
            ana_tab = rbind(ana_tab, record)
        }
    }
}
ana_tab
print(xtable(ana_tab), include.rownames=F)
