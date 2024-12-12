set.seed(2024)
library(rjags)

# number of subtrials
K <- 2

# observed data
n_k <- c(180,110)
N <- sum(n_k)
n_1k <- c(93,58)
basket <- rep(1:2,n_k)

covariate <- c("intolerant")
covariate <- c("haemoglobin")
covariate <- c("intolerant","haemoglobin")
covariate <- intersect(colnames(PV_centered),colnames(ET_centered))[1:10]
X_centered <- rbind(as.matrix(PV_centered[,covariate]),
                    as.matrix(ET_centered[,covariate]))
d <- dim(X_centered)[2]

T_ik <- c(PV$treatment,ET$treatment)
Y <- c(PV_CR,ET_CR)

# model
mod_file <- c("BHM.txt","adjBHM.txt","intBHM.txt")

# parameters and initial start
alpha <- rep(0,2)
beta <- matrix(rep(0,2*d),ncol=2)
gamma <- matrix(rep(0,2*d),ncol=2)

mu.theta <- 0
tau.theta <- 1
sigma.theta <- 1
re.theta <- rep(0,2)

theta <- mu.theta+re.theta

par <- list(BHM = c("alpha","theta","mu.theta","tau.theta"),
            adjBHM = c("alpha","beta","theta","mu.theta","tau.theta"),
            intBHM = c("alpha","beta","gamma","theta","mu.theta","tau.theta"))

ini <- list(BHM = list(alpha=alpha,re.theta=re.theta),
            adjBHM = list(alpha=alpha,beta=beta,re.theta=re.theta),
            intBHM = list(alpha=alpha,beta=beta,gamma=gamma,re.theta=re.theta))


# draw samples from the posterior distributions using MCMC
post <- function(K,N,basket,Y,T_ik,X,prec.HT,df.HT,mod_file,n.chains,ini,n.adapt,par,n.iter){
  if (mod_file=="BHM.txt") {
    data <- list(K=K,N=N,basket=basket,Y=Y,T_ik=T_ik,prec.HT=prec.HT,df.HT=df.HT)
  } else {
    data <- list(K=K,N=N,basket=basket,Y=Y,T_ik=T_ik,X=X,d=d,prec.HT=prec.HT,df.HT=df.HT)
  }
  
  mod <- jags.model(file=mod_file,data=data,n.chains = n.chains,inits = ini,n.adapt = n.adapt)
  post_sample <- coda.samples(mod,var=par,n.iter = n.iter)
  
  return(post_sample)
}

post_theta <- data.frame(matrix(ncol=3,nrow=10000))
colnames(post_theta) <- c("BHM","adjBHM","intBHM")
post_theta <- list(post_theta,post_theta)

theta_index <- sapply(1:K,function(k) paste0("theta[",k,"]"))

for (j in 1:3){
  post_sample <- post(K,N,basket,Y,T_ik,X_centered,prec.HT=1/2,df.HT=5,mod_file=mod_file[j],n.chains=1,ini=ini[[j]],n.adapt=2000,par=par[[j]],n.iter=10000)
  
  for (k in 1:K){
    post_theta[[k]][,j] <- post_sample[[1]][,theta_index[k]]
  }
}

###################################################################################
#####                                                                         #####                                                                       
#####                      Figure                                             #####
#####                                                                         ##### 
###################################################################################

library(reshape2)
library(ggplot2)
library(cowplot)

theta_index <- c()
theta_index[1]<-expression(theta[PV])
theta_index[2]<-expression(theta[ET])

figure<-list()


for (k in 1:2){
  post_theta_melt <- melt(post_theta[[k]])
  colnames(post_theta_melt) <- c("method","value")
  figure[[k]] <- ggplot(post_theta_melt)+
    stat_density(aes(x=value,color=method,linetype=method),position="identity",
                 geom="line",adjust = 1.5)+
    labs(title=theta_index[k])+
    xlab("")+
    theme(axis.text=element_text(size=10),
          axis.title.y = element_text(size=12),
          axis.text.y=element_text(size=10))
}

p <- plot_grid(figure[[1]],figure[[2]],ncol=2)
ggsave("postdensity_casestudy.pdf",plot=p,height = 4,width = 12)


###################################################################################
#####                                                                         #####                                                                       
#####                      Table                                              #####
#####                                                                         ##### 
###################################################################################

k <- 2
result_table <- data.frame(matrix(ncol=5,nrow=3))
colnames(result_table) <- c("mean","SD","CI lower","CI upper","SD reduction")
rownames(result_table) <- c("BHM","adjBHM","intBHM")

result_table[,1] <- colMeans(post_theta[[k]][,1:3])
result_table[,2] <- apply(post_theta[[k]][,1:3],2,sd)
result_table[,3] <- sapply(1:3, function(j) quantile(post_theta[[k]][,j],0.025))
result_table[,4] <- sapply(1:3, function(j) quantile(post_theta[[k]][,j],0.975))
result_table[,5] <- (1-result_table[,2]/result_table[1,2])*100

result_table <- round(result_table,3)
result_table[,5] <- paste0(result_table[,5],"%")

library(xtable)
xtable(result_table,digits=3)
