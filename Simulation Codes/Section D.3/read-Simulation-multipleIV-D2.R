library(ggplot2)
library(ggpubr)

bias_mat = rmse_mat = cover_mat = length_mat = qcomp_mat = NULL
for (a in seq(0.1, 1, 0.1)) {
    total = NULL
    for (round in seq(25)) {
        filename = paste0("results-pluralityRule/results",
                          "-a", a,
                          "-round", round, ".rds")
        one_round = readRDS(filename)
        total = rbind(total, one_round)
    }
    print(dim(total))
    beta0 = 1
    coef_tsht = coef_ciiv = coef_tsci = rep(NA, nrow(total))
    CI_tsht = CI_ciiv = CI_tsci = matrix(NA, nrow = nrow(total), ncol = 2)
    qcomp = rep(NA, nrow(total))
    for (i in 1:nrow(total)) {
        coef_tsht[i] = total[i,"res_tsht"][[1]]$betaHat[[1]]
        coef_ciiv[i] = total[i,"res_CIIV"][[1]]$Coefficients_CIM[1]
        if (class(total[i,"res_tsht"][[1]]$ci) == "list") {
            CI_tsht[i,] = total[i,"res_tsht"][[1]]$ci[[1]]
        } else {
            CI_tsht[i,] = total[i,"res_tsht"][[1]]$ci
        }
        CI_ciiv[i,] = total[i,"res_CIIV"][[1]]$ci_CIM
        qcomp[i] = total[i,"res_tsci"][[1]]$q.comp
        coef_tsci[i] = total[i,"res_tsci"][[1]]$Coef.vec[qcomp[i]+4]
        CI_tsci[i,] = c(
            coef_tsci[i] + qnorm(0.025)*total[i,"res_tsci"][[1]]$sd.vec[qcomp[i]+4],
            coef_tsci[i] + qnorm(0.975)*total[i,"res_tsci"][[1]]$sd.vec[qcomp[i]+4]
        )
    }
    bias_mat = rbind(bias_mat, c(a,
                                 mean(coef_tsht) - beta0,
                                 mean(coef_ciiv) - beta0,
                                 mean(coef_tsci) - beta0))
    rmse_mat = rbind(rmse_mat, c(a,
                                 sqrt(mean((coef_tsht - beta0)^2)),
                                 sqrt(mean((coef_ciiv - beta0)^2)),
                                 sqrt(mean((coef_tsci - beta0)^2))))
    cover_mat = rbind(cover_mat, c(a,
                                   mean(CI_tsht[,1] <= beta0 & beta0 <= CI_tsht[,2]),
                                   mean(CI_ciiv[,1] <= beta0 & beta0 <= CI_ciiv[,2]),
                                   mean(CI_tsci[,1] <= beta0 & beta0 <= CI_tsci[,2])))
    length_mat = rbind(length_mat, c(a,
                                     mean(CI_tsht[,2] - CI_tsht[,1]),
                                     mean(CI_ciiv[,2] - CI_ciiv[,1]),
                                     mean(CI_tsci[,2] - CI_tsci[,1])))
    qcomp_mat = rbind(qcomp_mat, c(a,
                                   sum(qcomp == 0),
                                   sum(qcomp == 1),
                                   sum(qcomp == 2)))
}
colnames(bias_mat) = colnames(rmse_mat) = colnames(cover_mat) = colnames(length_mat) = c("a","TSHT","CIIV","TSCI")
colnames(qcomp_mat) = c("a","0","1","2")

# plot --------------------------------------------------------------------

### RMSE
rmse_data = data.table(rmse_mat)
rmse_data = melt(rmse_data, id.vars=1)
colnames(rmse_data)[3] = "RMSE"
p_rmse = ggplot(data=rmse_data, aes(x=a, y=RMSE, col=variable)) +
    geom_line(lwd=1) +
    geom_point(size=1.8) +
    scale_x_continuous(breaks=seq(0.1, 1, 0.2)) +
    theme(
        legend.position = c(0.8,0.8)
        , legend.title = element_blank()
    ) +
    scale_color_manual(breaks=c("TSCI","TSHT","CIIV"),
                       values=c("#F8766D","#00BA38","#619CFF")) + 
    xlab("Invalidity level")

### coverage
cover_data = data.table(cover_mat)
cover_data = melt(cover_data, id.vars=1)
colnames(cover_data)[3] = "Coverage"
p_cover = ggplot(data=cover_data, aes(x=a, y=Coverage, col=variable)) +
    geom_line(lwd=1) +
    geom_point(size=1.8) +
    scale_x_continuous(breaks=seq(0.1, 1, 0.2)) +
    scale_y_continuous(breaks=seq(0,1,0.2), limits=c(0,1)) +
    geom_hline(yintercept = 0.95, 
               lwd = 0.5, 
               lty = 2) +
    theme(
        legend.position = c(0.8,0.3)
        , legend.title = element_blank()
    ) +
    scale_color_manual(breaks=c("TSCI","TSHT","CIIV"),
                       values=c("#F8766D","#00BA38","#619CFF")) + 
    xlab("Invalidity level")

### length
length_data = data.table(length_mat)
length_data = melt(length_data, id.vars=1)
colnames(length_data)[3] = "Length"
p_length = ggplot(data=length_data, aes(x=a, y=Length, col=variable)) +
    labs(col="") +
    geom_line(lwd=1) +
    geom_point(size=1.8) +
    scale_x_continuous(breaks=seq(0.1, 1, 0.2)) +
    scale_color_manual(breaks=c("TSCI","TSHT","CIIV"),
                       values=c("#F8766D","#00BA38","#619CFF")) + 
    theme(
        legend.position = c(0.8,0.8)
        , legend.title = element_blank()
    ) + 
    xlab("Invalidity level")

### q selection
qcomp_data = data.table(qcomp_mat)
qcomp_data = melt(qcomp_data, id.vars = 1)
qcomp_data$variable = factor(qcomp_data$variable, levels = c("2","1","0"))
p_qcomp = ggplot(qcomp_data, aes(fill = variable, y = value, x = a)) + 
    geom_bar(position = position_stack(), 
             width = 0.05,
             stat = "identity") +
    scale_x_continuous(breaks=seq(0.1, 1, 0.2)) +
    scale_fill_discrete(breaks = c("0","1","2"),
                        labels = c("q=0","q=1","q=2"),
                        type=c("#9E9E9E","#8983BF","#FFBE7A")) +
    theme(
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.key.size = unit(0.4, 'cm')
    ) +
    ylab("Count") + 
    xlab("Invalidity level")
    
p2 = ggarrange(ggarrange(p_rmse, p_cover, p_length, ncol=3, legend = "bottom", common.legend = T), 
               p_qcomp, 
               ncol=2, legend = "bottom", common.legend = F,
               widths = c(3, 1))
p2

ggarrange(p1, p2, nrow = 2, heights = c(0.8, 1), labels = c("D1","D2"), hjust = 0)
