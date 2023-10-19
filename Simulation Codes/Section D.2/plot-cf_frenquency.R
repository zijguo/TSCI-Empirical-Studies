# decomposition of TSLS point estimator: fhat as IV (with intercept)
source('Source-RF-hetero2.R')
library(ranger)
library(Matrix)
library(MASS)
library(ivreg)
library(ggplot2)
library(ggpubr)
library(data.table)

total = readRDS("results-MLIV_cf.rds")

# plot boxplot for coef_mliv and coef_tsci
tsci_init = unlist(total[,"tsci_nume"]) / unlist(total[,"tsci_deno"])
mliv_coef = unlist(total[,"mliv_nume"]) / unlist(total[,"mliv_deno"])
est_data = data.table(tsci = tsci_init, mliv = mliv_coef)
est_data = melt(est_data)
p_est = ggplot(data = est_data, aes(x = variable, y = value)) +
    geom_boxplot(width = 0.2) +
    # geom_violin() +
    coord_flip() +
    scale_x_discrete(labels = c('TSCI','MLIV')) +
    xlab("Methods") + 
    ylab("Point Estimators") + 
    theme_light()
p_est


# plot cf histograms
cf_data = data.table(unlist(total[,"c_f"]))
p_cf = ggplot(data = cf_data, aes(x = V1)) +
        geom_histogram(binwidth = 0.2, fill = "salmon", color = "black") + 
        xlab("c_f") + 
        ylab("Frequency") + 
        theme_light()
p_cf

ggarrange(p_est, p_cf, ncol = 2, widths = c(2,1))


