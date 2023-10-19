library(xtable)

ana_table = NULL
for (a in seq(0.15,0.35,0.02)) {
    total = NULL
    for (round in seq(25)) {
        filename = paste0("results-thol40/results",
                          "-a", a,
                          "-round", round, ".rds")
        one_round = readRDS(filename)
        total = rbind(total, one_round)
    }
    print(dim(total))
    beta0 = 0.5
    ivstr = mean(total[,3])
    rmse = sqrt(mean((total[,1] - beta0)^2)) 
    bias = mean(total[,1] - beta0)
    ci = cbind(total[,1] + qnorm(0.025)*total[,2],
               total[,1] + qnorm(0.975)*total[,2])
    cover = mean(ci[,1] <= beta0 & beta0 <= ci[,2])
    ana_table = rbind(ana_table, c(a,
                                   ivstr,
                                   rmse,
                                   bias,
                                   cover))
}
colnames(ana_table) = c("a","iv_str","rmse","bias","cover")


### coverage
cover_data = data.table(a = ana_table[order(ana_table[,"iv_str"]), "a"],
                        ivstr = ana_table[order(ana_table[,"iv_str"]), "iv_str"],
                        cover = ana_table[order(ana_table[,"iv_str"]), "cover"])
p_cover = ggplot(data = cover_data, aes(x = ivstr, y = cover)) + 
    geom_line(lwd = 0.7) + 
    geom_point() + 
    scale_x_continuous(breaks = seq(0,150,20)) + 
    scale_y_continuous(breaks = seq(0.3, 1, 0.1)) + 
    geom_hline(yintercept = 0.95, lty = 2) + 
    geom_vline(xintercept = 40, col = "red", lty = 4, lwd = 0.6) + 
    xlab("IV Strength") + 
    ylab("Coverage")
p_cover

### length
rmse_data = data.table(a = ana_table[order(ana_table[,"iv_str"]), "a"],
                         ivstr = ana_table[order(ana_table[,"iv_str"]), "iv_str"],
                         rmse = ana_table[order(ana_table[,"iv_str"]), "rmse"])
p_rmse = ggplot(data = rmse_data, aes(x = ivstr, y = rmse)) + 
    geom_line(lwd = 0.7) + 
    geom_point() + 
    scale_x_continuous(breaks = seq(0,150,20)) +
    geom_vline(xintercept = 40, col = "red", lty = 4, lwd = 0.6) + 
    xlab("IV Strength") + 
    ylab("RMSE")
p_rmse

### bias
bias_data = data.table(a = ana_table[order(ana_table[,"iv_str"]), "a"],
                         ivstr = ana_table[order(ana_table[,"iv_str"]), "iv_str"],
                         bias = ana_table[order(ana_table[,"iv_str"]), "bias"])
p_bias = ggplot(data = bias_data, aes(x = ivstr, y = bias)) + 
    geom_line(lwd = 0.7) + 
    geom_point() + 
    scale_x_continuous(breaks = seq(0,150,20)) +
    geom_vline(xintercept = 40, col = "red", lty = 4, lwd = 0.6) + 
    xlab("IV Strength") + 
    ylab("Bias")
p_bias

ggarrange(p_rmse, p_bias, p_cover, ncol = 3)

