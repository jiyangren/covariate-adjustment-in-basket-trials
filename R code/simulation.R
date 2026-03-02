library(rjags)
set.seed(2025)

# number of subtrials
K <- 4

# treatment effect
theta_k <- c(0,0,1,1)

# model
mod_file <- c("BHM.txt","adjBHM.txt","intBHM.txt")

# parameters and initial start
alpha <- rep(0,4)
beta <- rep(0,4)
gamma <- rep(0,4)
tau.y <- 1

mu.theta <- 0
tau.theta <- 1
sigma.theta <- 1
re.theta <- rep(0,4)

theta <- mu.theta+re.theta

par <- list(BHM = c("alpha","tau.y","theta","mu.theta","tau.theta"),
            adjBHM = c("alpha","beta","tau.y","theta","mu.theta","tau.theta"),
            intBHM = c("alpha","beta","gamma","tau.y","theta","mu.theta","tau.theta"))

ini <- list(BHM = list(alpha=alpha,tau.y=tau.y,re.theta=re.theta),
            adjBHM = list(alpha=alpha,beta=beta,tau.y=tau.y,re.theta=re.theta),
            intBHM = list(alpha=alpha,beta=beta,gamma=gamma,tau.y=tau.y,re.theta=re.theta))


# data generation process for scenario 1 2 3
scenario_1 <- function(K,n_k,N,n_1k,theta_k){
  X <- rnorm(N,0,1)
  X_centered <- c()
  T_ik <- rep(0,N)
  Y <- c()
  
  for (k in 1:K){
    if (k==1){
      T_ik[sample(1:n_k[1],n_1k[1])] <- 1
      X_centered[1:n_k[1]] <- X[1:n_k[1]]-mean(X[1:n_k[1]])
      Y[1:n_k[1]] <- 1+theta_k[k]*T_ik[1:n_k[1]]+4*(T_ik[1:n_k[1]]-1/2)*X_centered[1:n_k[1]]+rnorm(n_k[1],0,1)
    } else{
      unit <- (sum(n_k[1:(k-1)])+1):sum(n_k[1:k])
      T_ik[sample(unit,n_1k[k])] <- 1
      X_centered[unit] <- X[unit]-mean(X[unit])
      Y[unit] <- 1+theta_k[k]*T_ik[unit]+4*(T_ik[unit]-1/2)*X_centered[unit]+rnorm(n_k[k],0,1)
    }
  }
  return(list(Y=Y,T_ik=T_ik,X=X,X_centered=X_centered))
}

scenario_2 <- function(K,n_k,N,n_1k,theta_k){
  X <- rnorm(N,0,1)
  X_centered <- c()
  T_ik <- rep(0,N)
  Y <- c()
  
  for (k in 1:K){
    if (k==1){
      T_ik[sample(1:n_k[1],n_1k[1])] <- 1
      X_centered[1:n_k[1]] <- X[1:n_k[1]]-mean(X[1:n_k[1]])
      Y[1:n_k[1]] <- 1+theta_k[k]*T_ik[1:n_k[1]]-2*X_centered[1:n_k[1]]+rnorm(n_k[1],0,1)
    } else{
      unit <- (sum(n_k[1:(k-1)])+1):sum(n_k[1:k])
      T_ik[sample(unit,n_1k[k])] <- 1
      X_centered[unit] <- X[unit]-mean(X[unit])
      Y[unit] <- 1+theta_k[k]*T_ik[unit]-2*X_centered[unit]+rnorm(n_k[k],0,1)
    }
  }
  return(list(Y=Y,T_ik=T_ik,X=X,X_centered=X_centered))
}

scenario_3 <- function(K,n_k,N,n_1k,theta_k){
  X <- rnorm(N,0,1)
  X_centered <- c()
  T_ik <- rep(0,N)
  Y <- c()
  
  for (k in 1:K){
    if (k==1){
      T_ik[sample(1:n_k[1],n_1k[1])] <- 1
      X_centered[1:n_k[1]] <- X[1:n_k[1]]-mean(X[1:n_k[1]])
      Y[1:n_k[1]] <- 1+theta_k[k]*T_ik[1:n_k[1]]-2*X_centered[1:n_k[1]]+ T_ik[1:n_k[1]]*X_centered[1:n_k[1]]^3+rnorm(n_k[1],0,1)
    } else{
      unit <- (sum(n_k[1:(k-1)])+1):sum(n_k[1:k])
      T_ik[sample(unit,n_1k[k])] <- 1
      X_centered[unit] <- X[unit]-mean(X[unit])
      Y[unit] <- 1+theta_k[k]*T_ik[unit]-2*X_centered[unit]+T_ik[1:n_k[1]]*X_centered[unit]^3+rnorm(n_k[k],0,1)
    }
  }
  return(list(Y=Y,T_ik=T_ik,X=X,X_centered=X_centered))
}



# draw samples from the posterior distributions using MCMC
post <- function(K,N,basket,Y,T_ik,X,prec.HT,df.HT,mod_file,n.chains,ini,n.adapt,par,n.iter){
  if (mod_file=="BHM.txt") {
    data <- list(K=K,N=N,basket=basket,Y=Y,T_ik=T_ik,prec.HT=prec.HT,df.HT=df.HT)
  } else {
    data <- list(K=K,N=N,basket=basket,Y=Y,T_ik=T_ik,X=X,prec.HT=prec.HT,df.HT=df.HT)
  }

  mod <- jags.model(file=mod_file,data=data,n.chains = n.chains,inits = ini,n.adapt = n.adapt)
  post_sample <- coda.samples(mod,var=par,n.iter = n.iter)
  
  return(post_sample)
}


# calculate metrics from the posterior samples
metrics <- function(K,post_sample,n.adapt){
  result <- summary(window(post_sample,start=n.adapt))
  
  theta_index <- sapply(1:K,function(k) paste0("theta[",k,"]"))
  
  hat_theta <- c()
  rejection <- c()
  
  for (k in 1:K){
    hat_theta[k] <- result$statistics[theta_index[k],1]
    CI <- quantile(post_sample[[1]][,theta_index[k]],c(0.025,0.975))
    rejection[k] <- ifelse(CI[1]<0 & CI[2]>0,0,1)
  }
  
  return(list(hat_theta=hat_theta,rejection=rejection))
}


# simulate scenario 1 (or 2 or 3) 1000 times

cl <- makeCluster(4)

for (i in 10:1){
  
  # sample size, change from 56 to 560
  n_k <- c(8,12,16,20)*i
  N <- sum(n_k)
  
  # allocation ratio treatment:control = 1:1 or 3:1
  e <- 0.5 
  n_1k <-n_k*e
  
  # subtrial indicator
  basket <- rep(1:4,n_k)
  
  sim <- function(s){
    
    data <- scenario_1(K,n_k,N,n_1k,theta_k)
    T_ik <- data$T_ik
    Y <- data$Y
    X <- data$X
    X_centered <- data$X_centered
    
    hat_theta <- data.frame(matrix(ncol=12,nrow=4))
    colnames(hat_theta) <- c("BHM","adjBHM","intBHM","BHM2","adjBHM2","intBHM2","sANOVA","sANCOVAI","sANCOVAII","reANOVA","reANCOVAI","reANCOVAII")
    
    rejection <- hat_theta
    
    mu_theta <- c()
    power <- c()
    
    # BHMs
    for (j in 1:3){
      # draw posterior sample
      post_sample <- post(K,N,basket,Y,T_ik,X_centered,prec.HT=1/2,df.HT=5,mod_file=mod_file[j],n.chains=1,ini=ini[[j]],n.adapt=2000,par=par[[j]],n.iter=10000)
      # calculate metrics from the posterior sample
      metrics_result <- metrics(K,post_sample,n.adapt=2000)
      
      for (k in 1:K){
        hat_theta[k,j]<-metrics_result$hat_theta[k]
        rejection[k,j]<-metrics_result$rejection[k]
      }
      
      mu_theta[j]<-metrics_result$mu_theta
      power[j]<-metrics_result$power
    }
    
    for (j in 4:6){
      # draw posterior sample
      post_sample <- post(K,N,basket,Y,T_ik,X_centered,prec.HT=1/100,df.HT=5,mod_file=mod_file[j-3],n.chains=1,ini=ini[[j-3]],n.adapt=2000,par=par[[j-3]],n.iter=10000)
      # calculate metrics from the posterior sample
      metrics_result <- metrics(K,post_sample,n.adapt=2000)
      
      for (k in 1:K){
        hat_theta[k,j]<-metrics_result$hat_theta[k]
        rejection[k,j]<-metrics_result$rejection[k]
      }
      
      mu_theta[j]<-metrics_result$mu_theta
      power[j]<-metrics_result$power
    }
    
    # stratified ANCOVA
    sANOVA <- data.frame(basket=c(1:4),beta_T = c(1:4), se_T = c(1:4))
    sANCOVAI <- sANOVA
    sANCOVAII <- sANOVA
    
    for (k in 1:4){
      data_k <- data.frame(Y=Y[basket==k],T_ik=T_ik[basket==k],X_centered=X_centered[basket==k])
      
      model <- lm(Y ~ T_ik, data=data_k)
      sANOVA[k,2] <- coef(model)["T_ik"]
      sANOVA[k,3] <- sqrt(diag(vcov(model)))["T_ik"]
      
      model <- lm(Y ~ T_ik + X_centered, data=data_k)
      sANCOVAI[k,2] <- coef(model)["T_ik"]
      sANCOVAI[k,3] <- sqrt(diag(vcov(model)))["T_ik"]
      
      model <- lm(Y ~ T_ik*X_centered, data=data_k)
      sANCOVAII[k,2] <- coef(model)["T_ik"]
      sANCOVAII[k,3] <- sqrt(diag(vcov(model)))["T_ik"]
    }
    
    hat_theta[,7] <- sANOVA$beta_T
    hat_theta[,8] <- sANCOVAI$beta_T
    hat_theta[,9] <- sANCOVAII$beta_T
    
    rejection[,7] <- 1*(sANOVA$beta_T-qnorm(0.975)*sANOVA$se_T>0 | sANOVA$beta_T+qnorm(0.975)*sANOVA$se_T<0)
    rejection[,8] <- 1*(sANCOVAI$beta_T-qnorm(0.975)*sANCOVAI$se_T>0 | sANCOVAI$beta_T+qnorm(0.975)*sANCOVAI$se_T<0)
    rejection[,9] <- 1*(sANCOVAII$beta_T-qnorm(0.975)*sANCOVAII$se_T>0 | sANCOVAII$beta_T+qnorm(0.975)*sANCOVAII$se_T<0)
    
    sANOVA$weight <- 1 / (sANOVA$se_T^2)
    sANCOVAI$weight <- 1 / (sANCOVAI$se_T^2)
    sANCOVAII$weight <- 1 / (sANCOVAII$se_T^2)
    
    mu_theta[7] <- sum(sANOVA$beta_T * sANOVA$weight) / sum(sANOVA$weight)
    mu_theta[8] <- sum(sANCOVAI$beta_T * sANCOVAI$weight) / sum(sANCOVAI$weight)
    mu_theta[9] <- sum(sANCOVAII$beta_T * sANCOVAII$weight) / sum(sANCOVAII$weight)
    
    power[7] <- 1*(mu_theta[7]-qnorm(0.975)/sqrt(sum(sANOVA$weight))>0 | mu_theta[4]+qnorm(0.975)/sqrt(sum(sANOVA$weight))<0)
    power[8] <- 1*(mu_theta[8]-qnorm(0.975)/sqrt(sum(sANCOVAI$weight))>0 | mu_theta[5]+qnorm(0.975)/sqrt(sum(sANCOVAI$weight))<0)
    power[9] <- 1*(mu_theta[9]-qnorm(0.975)/sqrt(sum(sANCOVAII$weight))>0 | mu_theta[6]+qnorm(0.975)/sqrt(sum(sANCOVAII$weight))<0)
    
    
    # random effect ANCOVA
    reANOVA <- tryCatch(summary(lme(fixed = Y ~ T_ik,
                                    random = ~ T_ik | factor(basket),
                                    control = lmeControl(opt = "optim",maxIter=200))),
                        error = function(e) NA)
    reANCOVAI <- tryCatch(summary(lme(fixed = Y ~ T_ik + X_centered:factor(basket),
                                      random = ~ T_ik | factor(basket), 
                                      control = lmeControl(opt = "optim",maxIter=200))),
                          error = function(e) NA)
    reANCOVAII <- tryCatch(summary(lme(fixed = Y ~ T_ik * X_centered:factor(basket),
                                       random = ~ T_ik | factor(basket),
                                       control = lmeControl(opt = "optim",maxIter=200))),
                           error = function(e) NA)
    
    if (typeof(reANOVA)=="list"){
      mu_theta[10] <- reANOVA$tTable["T_ik",1]
      power[10] <- 1*(mu_theta[10]-qnorm(0.975)*reANOVA$tTable["T_ik",2]>0 | mu_theta[10]+qnorm(0.975)*reANOVA$tTable["T_ik",2]<0)
      hat_theta[,10] <- reANOVA$coefficients$random$`factor(basket)`[,"T_ik"]+mu_theta[10]
      var_re <- as.numeric(VarCorr(reANOVA)["T_ik","Variance"])
      sd_total <- sqrt(reANOVA$tTable["T_ik",2]^2 + var_re)
      rejection[,10] <- 1*(hat_theta[,10]-qnorm(0.975)*sd_total>0 | hat_theta[,10]+qnorm(0.975)*sd_total<0)
    } else {
      mu_theta[10]<-NA
      power[10]<-NA
      hat_theta[,10]<-NA
      rejection[,10]<-NA
    }
    
    if (typeof(reANCOVAI)=="list"){
      mu_theta[11] <- reANCOVAI$tTable["T_ik",1]
      power[11] <- 1*(mu_theta[11]-qnorm(0.975)*reANCOVAI$tTable["T_ik",2]>0 | mu_theta[11]+qnorm(0.975)*reANCOVAI$tTable["T_ik",2]<0)
      hat_theta[,11] <- reANCOVAI$coefficients$random$`factor(basket)`[,"T_ik"]+mu_theta[11]
      var_re <- as.numeric(VarCorr(reANCOVAI)["T_ik","Variance"])
      sd_total <- sqrt(reANCOVAI$tTable["T_ik",2]^2 + var_re)
      rejection[,11] <- 1*(hat_theta[,11]-qnorm(0.975)*sd_total>0 | hat_theta[,11]+qnorm(0.975)*sd_total<0)
    } else {
      mu_theta[11]<-NA
      power[11]<-NA
      hat_theta[,11]<-NA
      rejection[,11]<-NA
    }
    
    if (typeof(reANCOVAII)=="list"){
      mu_theta[12] <- reANCOVAII$tTable["T_ik",1]
      power[12] <- 1*(mu_theta[12]-qnorm(0.975)*reANCOVAII$tTable["T_ik",2]>0 | mu_theta[12]+qnorm(0.975)*reANCOVAII$tTable["T_ik",2]<0)
      hat_theta[,12] <- reANCOVAII$coefficients$random$`factor(basket)`[,"T_ik"]+mu_theta[12]
      var_re <- as.numeric(VarCorr(reANCOVAII)["T_ik","Variance"])
      sd_total <- sqrt(reANCOVAII$tTable["T_ik",2]^2 + var_re)
      rejection[,12] <- 1*(hat_theta[,12]-qnorm(0.975)*sd_total>0 | hat_theta[,12]+qnorm(0.975)*sd_total<0)
    } else {
      mu_theta[12]<-NA
      power[12]<-NA
      hat_theta[,12]<-NA
      rejection[,12]<-NA
    }
    
    return(list(hat_theta=hat_theta,rejection=rejection,
                mu_theta=mu_theta,power=power))
  }

  clusterExport(cl, varlist=c("K","n_k","N","n_1k","theta_k","basket","mod_file","ini","par",
                              "metrics","sim","scenario_1","scenario_2","scenario_3","post",
                              "jags.model","coda.samples",
                              "lme","lmeControl","VarCorr"))
  
  results <- parLapply(cl, 1:1000, sim)

  hat_theta <- data.frame(matrix(ncol=12,nrow=1000))
  colnames(hat_theta) <- c("BHM","adjBHM","intBHM","BHM2","adjBHM2","intBHM2","sANOVA","sANCOVAI","sANCOVAII","reANOVA","reANCOVAI","reANCOVAII")
  hat_theta <- list(hat_theta,hat_theta,hat_theta,hat_theta)
  
  rejection <- hat_theta
  
  for (k in 1:4){
    hat_theta[[k]] <- do.call(rbind,lapply(1:1000,function(s) results[[s]]$hat_theta[k,]))
    rejection[[k]] <- do.call(rbind,lapply(1:1000,function(s) results[[s]]$rejection[k,]))
  }
  
  save(results,hat_theta,rejection,file=paste0("result/scenario1_0.5/result_",i,".Rdata"))
}
