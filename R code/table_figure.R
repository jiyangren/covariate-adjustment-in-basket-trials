library(reshape2)
library(ggplot2)
library(cowplot)
library(xtable)

# run simulation.R first and store simulation data
# setwd("scenario1")


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


figure<-list()
for (k in 1:2){
  hat_theta_melt <- melt(hat_theta[[k]])
  colnames(hat_theta_melt) <- c("method","value")
  figure[[k]] <- ggplot(hat_theta_melt)+
    stat_density(aes(x=value,color=method,linetype=method),position="identity",
                 geom="line",adjust = 1.5)+
    geom_vline(aes(xintercept=0), linetype="dashed")+
    labs(title=theta_index[k])+
    xlab("")+
    theme(axis.text=element_text(size=10),
          axis.title.y = element_text(size=12),
          axis.text.y=element_text(size=10))
}

for (k in 3:4){
  hat_theta_melt<-melt(hat_theta[[k]])
  colnames(hat_theta_melt)<-c("method","value")
  figure[[k]] <- ggplot(hat_theta_melt)+
    stat_density(aes(x=value,color=method,linetype=method),position="identity",geom="line",adjust = 1.5)+
    geom_vline(aes(xintercept=1), linetype="dashed")+
    labs(title=theta_index[k])+
    xlab("")+
    theme(axis.text=element_text(size=10),
          axis.title.y = element_text(size=12),
          axis.text.y=element_text(size=10))
}
p<-plot_grid(figure[[1]],figure[[2]],figure[[3]],figure[[4]],ncol=2)

ggsave("density_scenario1.pdf",plot=p,height = 8,width = 12)


###################################################################################
#####                                                                         #####                                                                       
#####                      Table 3                                            #####
#####                                                                         ##### 
###################################################################################

k <- 1
result_table <- data.frame(matrix(ncol=5,nrow=3))
colnames(result_table) <- c("Bias","SD","RMSE","RMSE reduction","rejection rate")
rownames(result_table) <- c("BHM","adjBHM","intBHM")

result_table[,1] <- colMeans(hat_theta[[k]][,1:3])-theta_k[k]
result_table[,2] <- apply(hat_theta[[k]][,1:3],2,sd)
result_table[,3] <- sqrt(result_table[,1]^2+result_table[,2]^2)
result_table[,4] <- (1-result_table[,3]/result_table[1,3])*100
result_table[,5] <- colMeans(rejection[[k]][,1:3])

library(xtable)
xtable(result_table,digits=3)


###################################################################################
#####                                                                         #####                                                                       
#####                      Figure 4                                           #####
#####                                                                         ##### 
###################################################################################

theta_k<-c(0,0,1,1)

bias <- data.frame(matrix(ncol=4,nrow=10))
colnames(bias)<-c("size","BHM","adjBHM","intBHM")
bias[,1]<-c(1:10)*56
SD <- bias
RMSE <- bias
RMSE_reduction <- bias
typeIerror <- bias
power <- bias


for (i in 1:10){
  load(file=paste0("result_",i,".Rdata"))
  
  for (j in 1:3){
    bias[i,j+1] <- mean(sapply(1:4,function(k) abs(mean(hat_theta[[k]][,j])-theta_k[k])))
    SD[i,j+1] <- mean(sapply(1:4,function(k) sd(hat_theta[[k]][,j])))
    RMSE[i,j+1] <- mean(sapply(1:4,function(k) sqrt((mean(hat_theta[[k]][,j])-theta_k[k])^2+var(hat_theta[[k]][,j]))))
    RMSE_reduction[i,j+1] <- mean(sapply(1:4,function(k) 1-sqrt((mean(hat_theta[[k]][,j])-theta_k[k])^2+var(hat_theta[[k]][,j]))/sqrt((mean(hat_theta[[k]][,1])-theta_k[k])^2+var(hat_theta[[k]][,1]))))
    typeIerror[i,j+1] <- mean((rejection[[1]][,j]+rejection[[2]][,j])>=1)
    power[i,j+1] <- mean((rejection[[3]][,j]+rejection[[4]][,j])>=1)
  }
}

bias_melt <- melt(bias,id.vars="size")
bias_melt$metrics <- "bias"
SD_melt <- melt(SD,id.vars="size")
SD_melt$metrics <- "SD"
RMSE_melt <- melt(RMSE,id.vars="size")
RMSE_melt$metrics <- "RMSE"
RMSEreduction_melt <- melt(RMSE_reduction,id.vars="size")
RMSEreduction_melt$metrics <- "aversge RMSE reduction"
typeIerror_melt <- melt(typeIerror,id.vars="size")
typeIerror_melt$metrics <- "(family wise) type I error rate"
power_melt <- melt(power,id.vars="size")
power_melt$metrics <- "(disjunctive) power"

result <- rbind(RMSEreduction_melt,typeIerror_melt,power_melt)
result$metrics <- factor(result$metrics,levels=c("aversge RMSE reduction","(family wise) type I error rate","(disjunctive) power"))
colnames(result) <- c("size","method","value","metrics")
result$benchmark <- c()
result[which(result$metrics=="(family wise) type I error rate")[1],"benchmark"] <- 0.05
result[which(result$metrics=="(disjunctive) power")[1],"benchmark"] <- 0.8


p <- ggplot(data=result)+
  geom_hline(aes(yintercept=benchmark),linetype="dashed")+
  geom_text(aes(x=560,y=benchmark,label=benchmark))+
  geom_line(aes(x=size,y=value,color=method,linetype=method))+
  facet_wrap(vars(metrics),ncol=3)+
  xlab("sample size")+
  ylab("")+
  theme(axis.text=element_text(size=10),
        axis.text.y=element_text(size=10),
        strip.text=element_text(size=12))

ggsave("result_scenario1.pdf",plot=p,height = 4,width = 12)




