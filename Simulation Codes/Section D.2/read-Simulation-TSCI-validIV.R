library(xtable)
library(ggplot2)
library(data.table)
library(ggbreak)
library(ggpubr)
library(DescTools)  # for winsorized mean

# table -------------------------------------------------------------------

ana_tab_noself = NULL
for (a in seq(0.1,0.4,0.05)) {
    for (n in c(1000, 3000, 5000)) {
        total = NULL
        for (round in seq(25)) {
            filename = paste0("results-main/results-Simulation-TSCI-validIV-wtSelf0"
                              , "-a", a
                              , "-n", n
                              , "-round", round, ".rds")
            one_round = readRDS(filename)
            total = rbind(total, one_round)
        }
        print(dim(total))
        res_tsci = total[,"res_tsci"]
        res_mliv = total[,"res_mliv"]
        coef_tsci = coef_tsci0 = ivstr = coef_mliv = NULL
        for (i in 1:nrow(total)) {
            coef_tsci = c(coef_tsci, res_tsci[[i]][1])
            coef_tsci0 = c(coef_tsci0, res_tsci[[i]][2])
            ivstr = c(ivstr, res_tsci[[i]][3])
            coef_mliv = c(coef_mliv, res_mliv[[i]])
        }
        record = c(
            a = a,
            n = n,
            bias_tsci = mean((coef_tsci)) - 0.5,
            bias_tsci0 = mean((coef_tsci0)) - 0.5,
            bias_mliv = mean((coef_mliv)) - 0.5,
            mse_ratio = sqrt(mean(((coef_mliv - 0.5)^2))) / sqrt(mean(((coef_tsci - 0.5)^2))),
            ivstr = mean((ivstr))
        )
        ana_tab_noself = rbind(ana_tab_noself, record)
    }
}
ana_tab_wtself = NULL
for (a in seq(0.1,0.4,0.05)) {
    for (n in c(1000, 3000, 5000)) {
        total = NULL
        for (round in seq(25)) {
            filename = paste0("results-main/results-Simulation-TSCI-validIV-wtSelf1"
                              , "-a", a
                              , "-n", n
                              , "-round", round, ".rds")
            one_round = readRDS(filename)
            total = rbind(total, one_round)
        }
        print(dim(total))
        res_tsci = total[,"res_tsci"]
        res_mliv = total[,"res_mliv"]
        coef_tsci = coef_tsci0 = ivstr = coef_mliv = NULL
        for (i in 1:nrow(total)) {
            coef_tsci = c(coef_tsci, res_tsci[[i]][1])
            coef_tsci0 = c(coef_tsci0, res_tsci[[i]][2])
            ivstr = c(ivstr, res_tsci[[i]][3])
            coef_mliv = c(coef_mliv, res_mliv[[i]])
        }
        record = c(
            a = a,
            n = n,
            bias_tsci = mean((coef_tsci)) - 0.5,
            bias_tsci0 = mean((coef_tsci0)) - 0.5,
            bias_mliv = mean((coef_mliv)) - 0.5,
            mse_ratio = sqrt(mean(((coef_mliv - 0.5)^2))) / sqrt(mean(((coef_tsci - 0.5)^2))),
            ivstr = mean((ivstr))
        )
        ana_tab_wtself = rbind(ana_tab_wtself, record)
    }
}
ana_tab = cbind(ana_tab_noself, ana_tab_wtself[,-1:-2])
print(xtable(ana_tab, digits = c(0,2,0,rep(2,10))), include.rownames = F)

