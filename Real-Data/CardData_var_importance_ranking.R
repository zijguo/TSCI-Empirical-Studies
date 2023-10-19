library(ivmodel)
library(ggplot2)
source("Source-RF-hetero2.R", encoding="UTF-8")
data(card.data)

# random forests ----------------------------------------------------------

Xname2=c("exper", "expersq", "black", "south", "smsa", "reg661","reg662"
         , "reg663", "reg664", "reg665", "reg666", "reg667","reg668", "smsa66")

res <- NULL
for (i in 1:500) {
    print(i)
    ind = sample(1:nrow(card.data), nrow(card.data), replace = FALSE)
    Y = as.matrix(card.data[ind, "lwage"])
    D = as.matrix(card.data[ind, "educ"])
    Z = as.matrix(card.data[ind, "nearc4"])
    X = as.matrix(card.data[ind, Xname2])
    rf.result = TSCI.RF.fit(D, Z, X, mtry=1:(ncol(X)+1), min.node.size=c(3,5,10), num.trees=500, max.depth=0, split.prop=2/3, MSE.thol=1e12, forest.save=TRUE)  
    res = rbind(res, rf.result$forest$variable.importance)
}
colMeans(res)
saveRDS(res, "results-CardData_var_importance_ranking.rds")

# plot bar chart ----------------------------------------------------------

res = readRDS("results-CardData_var_importance_ranking.rds")
colnames(res) = c("nearc4", Xname2)
data_plot = data.frame(variable = colnames(res)
                       , importance = colMeans(res))
data_plot["class"] = c("IV", rep("Covariates",14))
data_plot = data_plot[order(data_plot$importance, decreasing=T),]

ggplot(data_plot, aes(x=reorder(variable, importance), y=importance, fill=class)) +
    geom_bar(stat='identity', width=0.7) +
    scale_fill_manual(values=c("steelblue", "salmon")) +
    coord_flip() +
    ylab("Variable Importance in Random Forests") +
    xlab(NULL) + 
    scale_y_continuous(breaks=seq(0,1800,200)) + 
    theme(legend.title=element_blank()
          , legend.position=c(0.8,0.2))
    
    


