library(rjags)
set.seed(2024)

# number of subtrials
K <- 4

# treatment effect
theta_k <- c(0,0,1,1)

# model
mod_file <- c("BHM","adjBHM","intBHM")

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
  T_ik <- c()
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
  T_ik <- c()
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
  T_ik <- c()
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

for (i in 1:10){
  
  # sample size, change from 56 to 560
  n_k <- c(8,12,16,20)*i
  N <- sum(n_k)
  
  # allocation ratio treatment:control = 1:1 or 3:1
  e <- 1/2 
  n_1k <-n_k*e
  
  # subtrial indicator
  basket <- rep(1:4,n_k)
  
  hat_theta <- data.frame(matrix(ncol=3,nrow=1000))
  colnames(hat_theta) <- c("BHM","adjBHM","intBHM")
  hat_theta <- list(hat_theta,hat_theta,hat_theta,hat_theta)
  rejection <- hat_theta
  
  # for each simulation loop
  for (s in 1:1000){
    # generate data
    data <- scenario_1(K,n_k,N,n_1k,theta_k)
    T_ik <- data$T_ik
    Y <- data$Y
    X <- data$X
    X_centered <- data$X_centered
    
    # for each model
    for (j in 1:3){
      # draw posterior sample
      post_sample <- post(K,N,basket,Y,T_ik,X_centered,prec.HT=1/2,df.HT=5,mod_file=mod_file[j],n.chains=1,ini=ini[[j]],n.adapt=2000,par=par[[j]],n.iter=10000)
      # calculate metrics from the posterior sample
      metrics_result <- metrics(K,post_sample,n.adapt=2000)
      
      for (k in 1:K){
        hat_theta[[k]][s,j]<-metrics_result$hat_theta[k]
        rejection[[k]][s,j]<-metrics_result$rejection[k]
      }
    }
  }
  
  save(hat_theta,rejection,file=paste0("scenario1/result_",i,".Rdata"))
}