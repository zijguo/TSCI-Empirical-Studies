# read DML regsDML setting
library(data.table)
library(ggplot2)
library(ggpubr)
library(DescTools)  # for winsorized mean
library(scales)  # for log scales

bias_mat = sd_mat = cover_mat = length_mat = rmse_mat = ivstr_mat = NULL
for (a in seq(0.6, 2, 0.2)) {
    total = NULL
    for (round in seq(25)) {
        filename = paste0("results-regsdmlSetting/results",
                          "-a", a,
                          "-round", round, ".rds")
        one_round = readRDS(filename)
        total = rbind(total, one_round)
    }
    print(dim(total))
    beta0 = 1
    ### bias
    bias_mat = rbind(bias_mat, c(a,
                                 # mean(unlist(total[,"beta_DML"]) - beta0),
                                 mean(unlist(total[,"beta_DMLorg"]) - beta0),
                                 # mean(unlist(total[,"beta_regsdml"]) - beta0),
                                 # mean(unlist(total[,"beta_tsls"]) - beta0),
                                 mean(unlist(total[,"beta_tsci"]) - beta0)
                                 ))
    ### coverage
    CI_DML = cbind(unlist(total[,"beta_DML"]) + qnorm(0.025)*unlist(total[,"sd_DML"]),
                   unlist(total[,"beta_DML"]) + qnorm(0.975)*unlist(total[,"sd_DML"]))
    CI_DMLorg = cbind(unlist(total[,"beta_DMLorg"]) + qnorm(0.025)*unlist(total[,"sd_DMLorg"]),
                      unlist(total[,"beta_DMLorg"]) + qnorm(0.975)*unlist(total[,"sd_DMLorg"]))
    CI_regsdml = cbind(unlist(total[,"beta_regsdml"]) + qnorm(0.025)*unlist(total[,"sd_regsdml"]),
                       unlist(total[,"beta_regsdml"]) + qnorm(0.975)*unlist(total[,"sd_regsdml"]))
    CI_tsls = cbind(unlist(total[,"beta_tsls"]) + qnorm(0.025)*unlist(total[,"sd_tsls"]),
                    unlist(total[,"beta_tsls"]) + qnorm(0.975)*unlist(total[,"sd_tsls"]))
    CI_tsci = cbind(unlist(total[,"beta_tsci"]) + qnorm(0.025)*unlist(total[,"sd_tsci"]),
                    unlist(total[,"beta_tsci"]) + qnorm(0.975)*unlist(total[,"sd_tsci"]))
    cover_mat = rbind(cover_mat, c(a,
                                   # mean(CI_DML[,1] <= beta0 & beta0 <= CI_DML[,2]),
                                   mean(CI_DMLorg[,1] <= beta0 & beta0 <= CI_DMLorg[,2]),
                                   # mean(CI_regsdml[,1] <= beta0 & beta0 <= CI_regsdml[,2]),
                                   # mean(CI_tsls[,1] <= beta0 & beta0 <= CI_tsls[,2]),
                                   mean(CI_tsci[,1] <= beta0 & beta0 <= CI_tsci[,2])
                                   ))
    ### length
    length_mat = rbind(length_mat, c(a,
                                   # mean(CI_DML[,2] - CI_DML[,1]),
                                   mean(CI_DMLorg[,2] - CI_DMLorg[,1]),
                                   # mean(CI_regsdml[,2] - CI_regsdml[,1]),
                                   # mean(CI_tsls[,2] - CI_tsls[,1]),
                                   mean(CI_tsci[,2] - CI_tsci[,1])
                                   ))
    ### RMSE
    rmse_mat = rbind(rmse_mat, c(a,
                               # sqrt(mean((unlist(total[,"beta_DML"]) - beta0)^2)),
                               sqrt(mean((unlist(total[,"beta_DMLorg"]) - beta0)^2)),
                               # sqrt(mean((unlist(total[,"beta_regsdml"]) - beta0)^2)),
                               # sqrt(mean((unlist(total[,"beta_tsls"]) - beta0)^2)),
                               sqrt(mean((unlist(total[,"beta_tsci"]) - beta0)^2))
                               ))
    ### ivstr
    ivstr_mat = rbind(ivstr_mat, c(a,
                                   mean(unlist(total[,"concen_DML"])),
                                   # mean(unlist(total[,"concen_tsls"])),
                                   mean(unlist(total[,"ivstr"]))
                                   ))
}
colnames(bias_mat)= colnames(cover_mat) = colnames(length_mat) = colnames(rmse_mat) = c("a","DML","TSCI")
colnames(ivstr_mat) = c("a","DML","TSCI")

# plot (no log) ---------------------------------------------------------------

### bias
bias_data = data.table(bias_mat)
bias_data = melt(bias_data, id.vars = 1)
bias_data[,3] = abs(bias_data[,3])
colnames(bias_data)[3] = "Bias"
p_bias = ggplot(data = bias_data, aes(x = a, y = Bias, col = variable)) +
    geom_line(lwd = 1) +
    geom_point(size = 1.8) +
    scale_x_continuous(breaks = seq(0.6,2,0.2)) +
    theme(legend.position = c(0.8,0.5),
          legend.title = element_blank()) +
    ylab("Absolute value of bias") +
    scale_color_manual(values = c("#619CFF","#F8766D"),
                       labels = c("DML","TSCI"))

### RMSE
rmse_data = data.table(rmse_mat)
rmse_data = melt(rmse_data, id.vars = 1)
colnames(rmse_data)[3] = "RMSE"
p_rmse = ggplot(data = rmse_data, aes(x = a, y = RMSE, col = variable)) +
    geom_line(lwd = 1) +
    geom_point(size = 1.8) +
    scale_x_continuous(breaks = seq(0.6,2,0.2)) +
    scale_y_continuous(breaks = c(0.05, 0.25, 0.45, 0.65)) + 
    theme(legend.position = c(0.8,0.4),
          legend.title = element_blank()) +
    scale_color_manual(values = c("#F8766D","#619CFF"),
                       breaks = c("TSCI","DML"))

### coverage
cover_data = data.table(cover_mat)
cover_data = melt(cover_data, id.vars = 1)
colnames(cover_data)[3] = "Coverage"
p_cover = ggplot(data = cover_data, aes(x = a, y = Coverage, col = variable)) +
    geom_line(lwd = 1) +
    geom_point(size = 1.8) +
    scale_x_continuous(breaks = seq(0.6,2,0.2)) +
    scale_y_continuous(breaks=seq(0.6,1,0.05), limits = c(0.9,1)) +
    geom_hline(yintercept = 0.95, 
               lwd = 0.5, 
               lty = 2) +
    theme(legend.position = c(0.3,0.2),
          legend.title = element_blank()) +
    scale_color_manual(values = c("#F8766D","#619CFF"),
                       breaks = c("TSCI","DML"))

### length
length_data = data.table(length_mat)
length_data = melt(length_data, id.vars = 1)
colnames(length_data)[3] = "Length"
p_length = ggplot(data = length_data, aes(x = a, y = Length, col = variable)) +
    labs(col = "") +
    geom_line(lwd = 1) +
    geom_point(size = 1.8) +
    scale_x_continuous(breaks = seq(0.6,2,0.2)) +
    scale_color_manual(values = c("#F8766D","#619CFF"),
                       breaks = c("TSCI","DML")) +
    theme(legend.position = c(0.8,0.8),
          legend.title = element_blank())

### iv strength
ivstr_data = data.table(ivstr_mat)
ivstr_data = melt(ivstr_data, id.vars = 1)
colnames(ivstr_data)[3] = "Value"
p_iv = ggplot(data = ivstr_data, aes(x = a, y = Value, col = variable)) +
    geom_line(lwd = 1) +
    geom_point(size = 1.8) +
    # scale_y_continuous(breaks=seq(0,1500,200)) +
    scale_x_continuous(breaks = seq(0.6,2,0.2)) +
    ylab("IV Strength") +
    scale_color_manual(values = c("#F8766D","#619CFF"),
                       breaks = c("TSCI","DML")) +
    theme(legend.position=c(0.8,0.2), 
          legend.title = element_blank())

# ggarrange(ggarrange(p_bias, p_sd, p_rmse, ncol = 3, legend = "none"), 
#           ggarrange(p_cover, p_length, p_iv, ncol = 3, common.legend = TRUE, legend = "bottom"),
#           nrow = 2)
p1 = ggarrange(p_rmse, p_cover, p_length, ncol = 3, common.legend = TRUE, legend = "none")
p1

