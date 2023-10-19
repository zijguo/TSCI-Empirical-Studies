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
            length_oracle = length_comp = length_robust = length_tsls = length_init = NULL
            length_plug_hetero = length_full_hetero = NULL
            qcomp = invalidity = NULL
            for (i in 1:nrow(total)) {
                CI_oracle = c(
                    total[i,"res_tsci"][[1]]$Coef.vec[5+vio]+qnorm(0.025)*total[i,"res_tsci"][[1]]$sd.vec[5+vio],
                    total[i,"res_tsci"][[1]]$Coef.vec[5+vio]+qnorm(0.975)*total[i,"res_tsci"][[1]]$sd.vec[5+vio]
                )
                if (total[i,"res_tsci"][[1]]$weak.iv == T) {
                    qcomp = c(qcomp, 0)
                    qrobust = 0
                    invalidity = c(invalidity, 0)
                } else {
                    qcomp = c(qcomp, total[i,"res_tsci"][[1]]$q.comp)
                    qrobust = total[i,"res_tsci"][[1]]$q.robust
                    invalidity = c(invalidity, total[i,"res_tsci"][[1]]$invalidity)
                }
                CI_comp = c(
                    total[i,"res_tsci"][[1]]$Coef.vec[5+qcomp[i]]+qnorm(0.025)*total[i,"res_tsci"][[1]]$sd.vec[5+qcomp[i]],
                    total[i,"res_tsci"][[1]]$Coef.vec[5+qcomp[i]]+qnorm(0.975)*total[i,"res_tsci"][[1]]$sd.vec[5+qcomp[i]]
                )
                CI_robust = c(
                    total[i,"res_tsci"][[1]]$Coef.vec[5+qrobust]+qnorm(0.025)*total[i,"res_tsci"][[1]]$sd.vec[5+qrobust],
                    total[i,"res_tsci"][[1]]$Coef.vec[5+qrobust]+qnorm(0.975)*total[i,"res_tsci"][[1]]$sd.vec[5+qrobust]
                )
                CI_tsls = total[i,"res_tsls"][[1]][2:3]
                CI_init = c(
                    total[i,"res_tsci"][[1]]$Coef.vec[vio+1]+qnorm(0.025)*total[i,"res_tsci"][[1]]$sd.vec[vio+1],
                    total[i,"res_tsci"][[1]]$Coef.vec[vio+1]+qnorm(0.975)*total[i,"res_tsci"][[1]]$sd.vec[vio+1]
                )
                CI_plug_hetero = c(
                    total[i,"res_plug_hetero"][[1]]$betaHat+qnorm(0.025)*total[i,"res_plug_hetero"][[1]]$sd
                    , total[i,"res_plug_hetero"][[1]]$betaHat+qnorm(0.975)*total[i,"res_plug_hetero"][[1]]$sd
                )
                CI_full_hetero = c(
                    total[i,"res_full_hetero"][[1]]$betaHat+qnorm(0.025)*total[i,"res_full_hetero"][[1]]$sd
                    , total[i,"res_full_hetero"][[1]]$betaHat+qnorm(0.975)*total[i,"res_full_hetero"][[1]]$sd
                )
                length_oracle = c(length_oracle, CI_oracle[2] - CI_oracle[1])
                length_comp = c(length_comp, CI_comp[2] - CI_comp[1])
                length_robust = c(length_robust, CI_robust[2] - CI_robust[1])
                length_tsls = c(length_tsls, CI_tsls[2] - CI_tsls[1])
                length_init = c(length_init, CI_init[2] - CI_init[1])
                length_plug_hetero = c(length_plug_hetero, CI_plug_hetero[2] - CI_plug_hetero[1])
                length_full_hetero = c(length_full_hetero, CI_full_hetero[2] - CI_full_hetero[1])
            }
            record = c(
                vio = vio
                , a = a
                , n = n
                , mean(length_oracle)
                , mean(length_comp)
                , mean(length_robust)
                , mean(invalidity)
                , sum(qcomp==0) / nrow(total)
                , sum(qcomp==1) / nrow(total)
                , sum(qcomp==2) / nrow(total)
                , sum(qcomp==3) / nrow(total)
                , mean(length_tsls)
                , mean(length_init)
                , mean(length_plug_hetero)
                , mean(length_full_hetero)
            )
            ana_tab = rbind(ana_tab, record)
        }
    }
}
ana_tab
print(xtable(ana_tab, digits = c(0,0,1,0,rep(2,12))), include.rownames=F)
