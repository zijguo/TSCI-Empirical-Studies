### we consider the bias, coverage, length, MSE
library(ggplot2)
library(ggpubr)

bias_mat = cover_mat = length_mat = mse_mat = ivstr_mat = qcomp_mat = NULL
for (a in c(seq(0,1.2,0.1))) {
    total = NULL
    for (round in seq(25)) {
        filename = paste0("results-sin-littleHidden/results"
                          , "-a", a
                          , "-round", round, ".rds")
        one_round = readRDS(filename)
        total = rbind(total, one_round)
    }
    print(dim(total))
    res_tsci = total[,"res_tsci"]
    res_tsciba_hetero = total[,"res_tsciba_hetero"]
    res_tsls = total[,"res_tsls"]
    res_dml = total[,"res_dml"]
    coef_tsci = coef_dml = coef_tsciba_hetero = coef_tsls = NULL
    CI_tsci = CI_dml = CI_tsciba_hetero = CI_tsls = NULL
    ivstr = concen = NULL
    qcomp = NULL
    for (i in 1:nrow(total)) {
        qcomp = c(qcomp, res_tsci[[i]]$q.comp)
        coef_tsci = c(coef_tsci, res_tsci[[i]]$Coef.vec[qcomp[i]+5])
        coef_dml = c(coef_dml, res_dml[[i]][1])
        if (res_tsciba_hetero[[i]]$run.OLS == T |
            res_tsciba_hetero[[i]]$weak.iv == T) {
            qcomp_ba = 0
        } else {
            qcomp_ba = res_tsciba_hetero[[i]]$q.comp
        }
        coef_tsciba_hetero = c(coef_tsciba_hetero, res_tsciba_hetero[[i]]$Coef.vec[qcomp_ba+1])
        coef_tsls = c(coef_tsls, res_tsls[[i]][1])
        CI_tsci = rbind(CI_tsci, c(
            coef_tsci[i]+qnorm(0.025)*res_tsci[[i]]$sd.vec[qcomp[i]+1],
            coef_tsci[i]+qnorm(0.975)*res_tsci[[i]]$sd.vec[qcomp[i]+1]))
        CI_dml = rbind(CI_dml, c(
            coef_dml[i]+qnorm(0.025)*res_dml[[i]][2],
            coef_dml[i]+qnorm(0.975)*res_dml[[i]][2]
        ))
        CI_tsciba_hetero = rbind(CI_tsciba_hetero, c(
            coef_tsciba_hetero[i]+qnorm(0.025)*res_tsciba_hetero[[i]]$sd.vec[qcomp_ba+1],
            coef_tsciba_hetero[i]+qnorm(0.975)*res_tsciba_hetero[[i]]$sd.vec[qcomp_ba+1]
        ))
        CI_tsls = rbind(CI_tsls, c(
            coef_tsls[i]+qnorm(0.025)*res_tsls[[i]][2],
            coef_tsls[i]+qnorm(0.975)*res_tsls[[i]][2]
        ))
        ivstr = c(ivstr, res_tsci[[i]]$iv.str[qcomp[i]+1])
        concen = c(concen, res_tsls[[i]][3])
    }    
    bias_mat = rbind(bias_mat, c(a, 
                                 mean(coef_tsls)-0.5,
                                 mean(coef_dml)-0.5,
                                 mean(coef_tsciba_hetero)-0.5,
                                 mean(coef_tsci)-0.5
    ))
    cover_mat = rbind(cover_mat, c(
        a,
        mean(CI_tsls[,1]<=0.5 & 0.5<=CI_tsls[,2]),
        mean(CI_dml[,1]<=0.5 & 0.5<=CI_dml[,2]),
        mean(CI_tsciba_hetero[,1]<=0.5 & 0.5<=CI_tsciba_hetero[,2]),
        mean(CI_tsci[,1]<=0.5 & 0.5<=CI_tsci[,2])
    ))
    length_mat = rbind(length_mat, c(
        a,
        mean(CI_tsls[,2] - CI_tsls[,1]),
        mean(CI_dml[,2] - CI_dml[,1]),
        mean(CI_tsciba_hetero[,2] - CI_tsciba_hetero[,1]),
        mean(CI_tsci[,2] - CI_tsci[,1])
    ))
    mse_mat = rbind(mse_mat, c(
        a,
        sqrt(mean((coef_tsls-0.5)^2)),
        sqrt(mean((coef_dml-0.5)^2)),
        sqrt(mean((coef_tsciba_hetero-0.5)^2)),
        sqrt(mean((coef_tsci-0.5)^2))
    ))
    ivstr_mat = rbind(ivstr_mat, c(
        a,
        mean(ivstr),
        mean(concen)
    ))
    qcomp_mat = rbind(qcomp_mat, c(
        a,
        sum(qcomp == 0),
        sum(qcomp == 1),
        sum(qcomp == 2),
        sum(qcomp == 3)
    ))
}
colnames(bias_mat) = colnames(cover_mat) = colnames(length_mat) = colnames(mse_mat) = c("a","TSLS","DML","TSCI-ba","TSCI")
colnames(qcomp_mat) = c("a","0","1","2","3")


# plot --------------------------------------------------------------------

### bias
bias_data = data.table(bias_mat[,-4])
bias_data = melt(bias_data, id.vars=1)
bias_data[,3] = abs(bias_data[,3])
colnames(bias_data)[3] = "Bias"
p_bias = ggplot(data=bias_data, aes(x=a, y=Bias, col=variable)) +
    geom_line(lwd=1) +
    geom_point(size=1.8) +
    scale_x_continuous(breaks=seq(0,1.2,0.2)) +
    # scale_y_continuous(breaks=seq(0,0.2,0.05)) +
    theme(
        legend.position = c(0.8,0.8)
        , legend.title = element_blank()
    ) +
    ylab("Absolute value of bias") +
    scale_color_manual(values = c("#F8766D","#00BA38","#619CFF"),
                       breaks = c("TSCI","TSLS","DML"))

### coverage
cover_data = data.table(cover_mat[,-4])
cover_data = melt(cover_data, id.vars=1)
colnames(cover_data)[3] = "Coverage"
p_cover = ggplot(data=cover_data, aes(x=a, y=Coverage, col=variable)) +
    geom_line(data = cover_data[cover_data$variable=="TSLS",], lwd = 1.8) +
    geom_line(lwd=1) +
    geom_point(size=1.8) +
    scale_x_continuous(breaks=seq(0,1.2,0.2)) +
    # scale_y_continuous(breaks=seq(0,1,0.1), limits=c(0,1)) +
    geom_hline(yintercept = 0.95, 
               lwd = 0.5, 
               lty = 2) +
    theme(
        legend.position = c(0.8,0.3)
        , legend.title = element_blank()
    ) +
    scale_color_manual(values = c("#F8766D","#619CFF","#00BA38"),
                       breaks = c("TSCI","DML","TSLS"))

### length
length_data = data.table(length_mat[,-4])
length_data = melt(length_data, id.vars=1)
colnames(length_data)[3] = "Length"
p_length = ggplot(data=length_data, aes(x=a, y=Length, col=variable)) +
    labs(col="") +
    geom_line(lwd=1) +
    geom_point(size=1.8) +
    scale_x_continuous(breaks=seq(0,1.2,0.2)) +
    scale_y_continuous(breaks=seq(0,0.3,0.1)) +
    scale_color_manual(values = c("#F8766D","#619CFF","#00BA38"),
                       breaks = c("TSCI","DML","TSLS")) +
    theme(
        legend.position = c(0.8,0.8)
        , legend.title = element_blank()
    )

### MSE
mse_data = data.table(mse_mat[,-4])
mse_data = melt(mse_data, id.vars=1)
colnames(mse_data)[3] = "RMSE"
p_mse = ggplot(data=mse_data, aes(x=a, y=RMSE, col=variable)) +
    geom_line(lwd=1) +
    geom_point(size=1.8) +
    scale_x_continuous(breaks=seq(0,1.2,0.2)) +
    theme(
        legend.position = c(0.8,0.8)
        , legend.title = element_blank()
    ) +
    scale_color_manual(values = c("#F8766D","#619CFF","#00BA38"),
                       breaks = c("TSCI","DML","TSLS"))

### q selection
qcomp_data = data.table(qcomp_mat)
qcomp_data = melt(qcomp_data, id.vars = 1)
qcomp_data$variable = factor(qcomp_data$variable, levels = c("3","2","1","0"))
p_qcomp = ggplot(qcomp_data, aes(fill = variable, y = value, x = a)) + 
    geom_bar(position = position_stack(), 
             width = 0.05,
             stat = "identity") +
    scale_x_continuous(breaks=seq(0,1.2,0.2)) +
    scale_fill_discrete(type=c("#8ECFC9","#9E9E9E","#8983BF","#FFBE7A"),
                        breaks=c("0","1","2","3"),
                        labels=c("q=0","q=1","q=2","q=3")) +
    theme(
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.key.size = unit(0.4, 'cm')
    ) +
    ylab("Count")

p3 = ggarrange(ggarrange(p_mse, p_cover, p_length, ncol=3, legend = "bottom", common.legend = T), 
          p_qcomp, 
          ncol=2, legend = "bottom", common.legend = F,
          widths = c(3, 1))
p3

ggarrange(p1, p2, p3,
          heights = c(0.75,0.75,1),
          nrow = 3,
          ncol = 1,
          labels = c("C1","C2","C3"),
          hjust = 0,
          vjust = 1)
