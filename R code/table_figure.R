library(reshape2)
library(ggplot2)
library(cowplot)
library(xtable)

# run simulation.R first and store simulation data
# setwd("scenario1_0.5")


###################################################################################
#####                                                                         #####                                                                       
#####                      Figure 3                                           #####
#####                                                                         ##### 
###################################################################################

load(file="result_10.Rdata")

theta_index <- c()
theta_index[1]<-expression(theta[1])
theta_index[2]<-expression(theta[2])
theta_index[3]<-expression(theta[3])
theta_index[4]<-expression(theta[4])


model_colors <- c("BHMs" = "#F8766D",
                  "stratified model" = "#00BA38",
                  "random effect model" = "#619CFF")

adjust_linetypes <- c("unadjusted" = "solid",
                      "ANCOVAI" = "dashed",
                      "ANCOVAII" = "dotdash") 

figure<-list()
for (k in 1:2){
  hat_theta_melt <- melt(hat_theta[[k]][,c(1:3,8:9,11:12)],id.vars = NULL)
  colnames(hat_theta_melt) <- c("method","value")
  
  hat_theta_melt$adjustment <- ifelse(hat_theta_melt$method=="BHM","unadjusted",
                                 ifelse(hat_theta_melt$method %in% c("adjBHM","sANCOVAI","reANCOVAI"),"ANCOVAI","ANCOVAII"))
  hat_theta_melt$adjustment <- factor(hat_theta_melt$adjustment,levels=c("unadjusted","ANCOVAI","ANCOVAII"))
  
  hat_theta_melt$model <- ifelse(hat_theta_melt$method %in% c("BHM","adjBHM","intBHM"),"BHMs",
                                 ifelse(hat_theta_melt$method %in% c("sANCOVAI","sANCOVAII"),"stratified model","random effect model"))
  hat_theta_melt$model <- factor(hat_theta_melt$model,levels=c("BHMs","stratified model","random effect model"))
  
  
  method_cols <- setNames(model_colors[as.character(hat_theta_melt$model)], hat_theta_melt$method)
  method_ltys <- setNames(adjust_linetypes[as.character(hat_theta_melt$adjustment)], hat_theta_melt$method)
  
  
  figure[[k]] <- ggplot(hat_theta_melt)+
    stat_density(aes(x=value,color=method,linetype=method),position="identity",
                 geom="line",adjust = 1.5,show.legend = FALSE)+
    geom_vline(aes(xintercept=0), linetype="dashed",show.legend = FALSE)+
    scale_color_manual(name = "method", values = method_cols, breaks = names(method_cols)) +
    scale_linetype_manual(name = "method", values = method_ltys, breaks = names(method_ltys)) +
    labs(title=theta_index[k])+
    xlab("")+
    theme(axis.text=element_text(size=10),
          axis.title.y = element_text(size=12),
          axis.text.y=element_text(size=10))
}

for (k in 3:4){
  hat_theta_melt <- melt(hat_theta[[k]][,c(1:3,8:9,11:12)],id.vars = NULL)
  colnames(hat_theta_melt) <- c("method","value")
  hat_theta_melt$adjustment <- ifelse(hat_theta_melt$method=="BHM","unadjusted",
                                 ifelse(hat_theta_melt$method %in% c("adjBHM","sANCOVAI","reANCOVAI"),"ANCOVAI","ANCOVAII"))
  hat_theta_melt$adjustment <- factor(hat_theta_melt$adjustment,levels=c("unadjusted","ANCOVAI","ANCOVAII"))
  
  hat_theta_melt$model <- ifelse(hat_theta_melt$method %in% c("BHM","adjBHM","intBHM"),"BHMs",
                                 ifelse(hat_theta_melt$method %in% c("sANCOVAI","sANCOVAII"),"stratified model","random effect model"))
  hat_theta_melt$model <- factor(hat_theta_melt$model,levels=c("BHMs","stratified model","random effect model"))
  
  method_cols <- setNames(model_colors[as.character(hat_theta_melt$model)], hat_theta_melt$method)
  method_ltys <- setNames(adjust_linetypes[as.character(hat_theta_melt$adjustment)], hat_theta_melt$method)
  
  figure[[k]] <- ggplot(hat_theta_melt)+
    stat_density(aes(x=value,color=method,linetype=method),position="identity",
                geom="line",adjust = 1.5,show.legend = FALSE)+
    geom_vline(aes(xintercept=1), linetype="dashed",show.legend = FALSE)+
    scale_color_manual(name = "method", values = method_cols, breaks = names(method_cols)) +
    scale_linetype_manual(name = "method", values = method_ltys, breaks = names(method_ltys)) +
    labs(title=theta_index[k])+
    xlab("")+
    theme(axis.text=element_text(size=10),
          axis.title.y = element_text(size=12),
          axis.text.y=element_text(size=10))
}

p<-plot_grid(figure[[1]],figure[[2]],figure[[3]],figure[[4]],ncol=2)

ggsave("density_scenario1_0.5.pdf",plot=p,height = 8,width = 12)

p <- ggplot(hat_theta_melt)+
  stat_density(aes(x=value,color=method,linetype=method),position="identity",
               geom="line",adjust = 1.5)+
  geom_vline(aes(xintercept=0), linetype="dashed",show.legend = FALSE)+
  scale_color_manual(name = "method", values = method_cols, breaks = names(method_cols)) +
  scale_linetype_manual(name = "method", values = method_ltys, breaks = names(method_ltys)) +
  xlab("")+
  theme(axis.text=element_text(size=10),
        axis.title.y = element_text(size=12),
        axis.text.y=element_text(size=10),
        legend.position = "bottom")

# only show legend
leg <- ggdraw() + draw_plot(get_legend(p))

ggsave("result/scenario1_0.5/legend.pdf",plot=leg,height = 1,width = 12)

###################################################################################
#####                                                                         #####                                                                       
#####                      Table 3                                            #####
#####                                                                         ##### 
###################################################################################

k <- 1
result_table <- data.frame(matrix(ncol=5,nrow=12))
colnames(result_table) <- c("Bias","SD","RMSE","RMSE reduction","rejection rate")
rownames(result_table) <- c("BHM","adjBHM","intBHM","BHM2","adjBHM2","intBHM2","sANOVA","sANCOVAI","sANCOVAII","reANOVA","reANCOVAI","reANCOVAII")

theta_k <- c(0,0,1,1)
result_table[,1] <- colMeans(hat_theta[[k]][,1:12])-theta_k[k]
result_table[,2] <- apply(hat_theta[[k]][,1:12],2,sd)
theta_square <- (hat_theta[[k]][,1:12]-theta_k[k])^2
result_table[,3] <- sqrt(colMeans(theta_square))
result_table[,4] <- (1-result_table[,3]/result_table[1,3])*100
result_table[,5] <- colMeans(rejection[[k]][,1:12])

MCSE <- data.frame(matrix(ncol=5,nrow=12))
colnames(MCSE) <- c("Bias","SD","RMSE","RMSE reduction","rejection rate")
rownames(MCSE) <- c("BHM","adjBHM","intBHM","BHM2","adjBHM2","intBHM2","sANOVA","sANCOVAI","sANCOVAII","reANOVA","reANCOVAI","reANCOVAII")

MCSE[,1] <- apply(hat_theta[[k]],2,sd)/sqrt(1000)
MCSE[,2] <- apply(hat_theta[[k]],2,sd)/sqrt(2*(1000-1))
MCSE[,3] <- 1/2/result_table[,3]*apply(theta_square,2,sd)/sqrt(1000)
MCSE[1,4] <- 0
var_MSEAvsB <- apply(theta_square,2,var)/1000/mean(theta_square[,1])^2+
  colMeans(theta_square)^2*var(theta_square[,1])/1000/mean(theta_square[,1])^4-
  2*colMeans(theta_square)*cov(theta_square)[1,]/1000/mean(theta_square[,1])^3
MCSE[2:12,4] <- sqrt(var_MSEAvsB[2:12])*result_table[1,3]*50/result_table[2:12,3]
MCSE[,5] <- sqrt(result_table[,5]*(1-result_table[,5])/1000)

combine_table <- result_table
for (i in 1:12){
  for (j in 1:5){
    combine_table[i,j] <- paste0(round(result_table[i,j],3)," (",round(MCSE[i,j],3),")")
  }
}

library(xtable)
xtable(combine_table)


###################################################################################
#####                                                                         #####                                                                       
#####                      Figure 4                                           #####
#####                                                                         ##### 
###################################################################################

theta_k<-c(0,0,1,1)

bias <- data.frame(matrix(ncol=13,nrow=10))
colnames(bias)<-c("size","BHM","adjBHM","intBHM","BHM2","adjBHM2","intBHM2","sANOVA","sANCOVAI","sANCOVAII","reANOVA","reANCOVAI","reANCOVAII")
bias[,1]<-c(1:10)*56
SD <- bias
RMSE <- bias
RMSE_reduction <- bias
typeIerror <- bias
Power <- bias
converge <- bias

for (i in 1:10){
  load(file=paste0("result/scenario1_0.5/result_",i,".Rdata"))
  
  converge[i,2:13] <- colMeans(is.na(hat_theta[[1]]))
  
  for (j in 1:12){
    bias[i,j+1] <- mean(sapply(1:4,function(k) abs(mean(hat_theta[[k]][,j],na.rm=T)-theta_k[k])))
    SD[i,j+1] <- mean(sapply(1:4,function(k) sd(hat_theta[[k]][,j],na.rm=T)))
    RMSE[i,j+1] <- mean(sapply(1:4,function(k) sqrt(mean((hat_theta[[k]][,j]-theta_k[k])^2,na.rm=T))))
    RMSE_reduction[i,j+1] <- mean(sapply(1:4,function(k) 1-sqrt(mean((hat_theta[[k]][,j]-theta_k[k])^2,na.rm=T))/sqrt(mean((hat_theta[[k]][,1]-theta_k[k])^2,na.rm=T))))
    typeIerror[i,j+1] <- mean((rejection[[1]][,j]+rejection[[2]][,j])>=1,na.rm=T)
    Power[i,j+1] <- mean((rejection[[3]][,j]+rejection[[4]][,j])>=1,na.rm=T)
  }
}
converge_melt <- melt(converge,id.vars="size")
converge_melt$metrics <- "probability of non-convergence"
bias_melt <- melt(bias,id.vars="size")
bias_melt$metrics <- "bias"
SD_melt <- melt(SD,id.vars="size")
SD_melt$metrics <- "SD"
RMSE_melt <- melt(RMSE,id.vars="size")
RMSE_melt$metrics <- "RMSE"
RMSEreduction_melt <- melt(RMSE_reduction,id.vars="size")
RMSEreduction_melt$metrics <- "average RMSE reduction"
typeIerror_melt <- melt(typeIerror,id.vars="size")
typeIerror_melt$metrics <- "(family wise) type I error rate"
power_melt <- melt(Power,id.vars="size")
power_melt$metrics <- "(disjunctive) power"

result <- rbind(converge_melt,RMSEreduction_melt,typeIerror_melt,power_melt)
result$metrics <- factor(result$metrics,levels=c("probability of non-convergence","average RMSE reduction","(family wise) type I error rate","(disjunctive) power"))
colnames(result) <- c("size","method","value","metrics")
result$benchmark <- c()
result[which(result$metrics=="(family wise) type I error rate")[1],"benchmark"] <- 1-(1-0.05)^2
result[which(result$metrics=="(disjunctive) power")[1],"benchmark"] <- 0.8

result$adjustment <- ifelse(result$method %in% c("BHM","sANOVA","reANOVA"),"unadjusted",
                            ifelse(result$method %in% c("adjBHM","sANCOVAI","reANCOVAI"),"ANCOVAI","ANCOVAII"))
result$adjustment <- factor(result$adjustment,levels=c("unadjusted","ANCOVAI","ANCOVAII"))

result$model <- ifelse(result$method %in% c("BHM","adjBHM","intBHM"),"BHMs",
                       ifelse(result$method %in% c("sANOVA","sANCOVAI","sANCOVAII"),"stratified model","random effect model"))

result$model <- factor(result$model,levels=c("BHMs","stratified model","random effect model"))

method_cols <- setNames(model_colors[as.character(result$model)], result$method)
method_ltys <- setNames(adjust_linetypes[as.character(result$adjustment)], result$method)

p1 <- ggplot(data=result[result$method %in% c("BHM","adjBHM","intBHM","sANCOVAI","sANCOVAII","reANCOVAI","reANCOVAII") & 
                           result$metrics=="probability of non-convergence",])+
  geom_hline(aes(yintercept=benchmark),linetype="dashed",show.legend=FALSE)+
  geom_line(aes(x=size,y=value,color=method,linetype=method))+
  xlab("sample size")+
  ylab("probability of non-convergence")+
  scale_color_manual(name = "method", values = method_cols, breaks = names(method_cols),guide = "none") +
  scale_linetype_manual(name = "method", values = method_ltys, breaks = names(method_ltys),guide = "none") +
  theme(axis.text=element_text(size=10),
        axis.title.y = element_text(size=12),
        axis.text.y=element_text(size=10))

p2 <- ggplot(data=result[result$method %in% c("BHM","adjBHM","intBHM","sANCOVAI","sANCOVAII","reANCOVAI","reANCOVAII") & 
                           result$metrics=="average RMSE reduction",])+
  geom_hline(aes(yintercept=benchmark),linetype="dashed",show.legend=FALSE)+
  geom_line(aes(x=size,y=value,color=method,linetype=method))+
  xlab("sample size")+
  ylab("average RMSE reduction")+
  scale_color_manual(name = "method", values = method_cols, breaks = names(method_cols),guide = "none") +
  scale_linetype_manual(name = "method", values = method_ltys, breaks = names(method_ltys),guide = "none") +
  theme(axis.text=element_text(size=10),
        axis.title.y = element_text(size=12),
        axis.text.y=element_text(size=10))

p3 <- ggplot(data=result[result$method %in% c("BHM","adjBHM","intBHM","sANCOVAI","sANCOVAII","reANCOVAI","reANCOVAII") & 
                           result$metrics=="(family wise) type I error rate",])+
  geom_hline(aes(yintercept=benchmark),linetype="dashed",show.legend=FALSE)+
  geom_line(aes(x=size,y=value,color=method,linetype=method))+
  xlab("sample size")+
  ylab("(family wise) type I error rate")+
  scale_color_manual(name = "method", values = method_cols, breaks = names(method_cols),guide = "none") +
  scale_linetype_manual(name = "method", values = method_ltys, breaks = names(method_ltys),guide = "none") +
  theme(axis.text=element_text(size=10),
        axis.title.y = element_text(size=12),
        axis.text.y=element_text(size=10))

p4 <- ggplot(data=result[result$method %in% c("BHM","adjBHM","intBHM","sANCOVAI","sANCOVAII","reANCOVAI","reANCOVAII") & 
                           result$metrics=="(disjunctive) power",])+
  geom_hline(aes(yintercept=benchmark),linetype="dashed",show.legend=FALSE)+
  geom_line(aes(x=size,y=value,color=method,linetype=method))+
  xlab("sample size")+
  ylab("(disjunctive) power")+
  scale_color_manual(name = "method", values = method_cols, breaks = names(method_cols),guide = "none") +
  scale_linetype_manual(name = "method", values = method_ltys, breaks = names(method_ltys),guide = "none") +
  theme(axis.text=element_text(size=10),
        axis.title.y = element_text(size=12),
        axis.text.y=element_text(size=10))

p<-plot_grid(p1,p2,p3,p4,ncol=2)

ggsave("result/scenario1_0.5/result_scenario1_0.5.pdf",plot=p,height = 6,width = 12)





