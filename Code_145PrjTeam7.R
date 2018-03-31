setwd("~/Dropbox/zijide 145Prj")
library(BayesLogit)
library(MASS)
library(MCMCpack)

diabetes = read.table('Diabetes.txt', header = TRUE)[,-1]
#View(diabetes)
#multicollinearity between independent variables
pairs(data = diabetes[-ncol(diabetes)], age ~.,lower.panel=NULL)
cor(diabetes$hip,diabetes$waist)
cor(diabetes$weight,diabetes$waist) 
cor(diabetes$weight,diabetes$hip) 
cor(diabetes$BloodPressureHigh,diabetes$BloodPressureLow)
#average waist, hip and weight 
diabetes$weight = (diabetes$hip+diabetes$weight+diabetes$waist)/3
pairs(data = diabetes[-which(names(diabetes)%in%c("waist","hip",'diabetes'))], age ~.,lower.panel=NULL)
diabetes = diabetes[-which(names(diabetes)%in%c("waist","hip"))]
diabetes$location = as.integer(diabetes$location) - 1
diabetes$gender = as.integer(diabetes$gender) - 1
#######################
#gibbs_sampler function
#######################

gibbs_sampler = function(n_sample, burnin, data = diabetes, mi=1, C, b0){
  X = cbind(1,as.matrix(diabetes[-ncol(diabetes)]))
  Y = as.matrix(diabetes[ncol(diabetes)])
  #View(X)
  #View(Y)
  n_data = dim(diabetes)[1]
  n_beta = dim(diabetes)[2] #beta -1(itself) +1 coefficient 
  MLE_logit_model = glm(diabetes ~.,family = binomial(link='logit'), data = diabetes)
  B = C*diag(1,nrow= n_beta,ncol=n_beta)  #prior covariance matrix 
  #C large, prior has less information; when C is small, more information
  intial_beta = matrix(rep(0, n_beta),nrow = n_beta)
  m_i = 1 
  BETA = matrix(NA,nrow=n_sample,ncol=n_beta)
  BETA[1,] = intial_beta
  for (i in 2:n_sample){
    # Draw w's from Polya-Gamma dist.
    yy = abs(X%*%matrix(BETA[(i-1),],nrow=n_beta))
    set.seed(23123441)
    W= rpg(num = n_data, h = 1, z=c(yy)) #c has length n 
    SIGMA_w = solve(t(X) %*% diag(W) %*% X + solve(B))
    m_w = SIGMA_w %*% (t(X) %*% (Y - 0.5*matrix(rep(m_i, n_data),nrow = n_data)) + solve(B)%*%b0)  
    # Draw beta from p-variate Normal dist.
    BETA[i,] = mvrnorm(1, mu=m_w, Sigma=SIGMA_w)
  }#for 
  #return(BETA)
  
  "sum.stats" <- function(x){
    nstats <- 5
    ret <- rep(NA,nstats)
    names(ret) <- c("mean","sd","median","q02.5","q97.5")
    ret[1:2] <- c(mean(x),sd(x))
    ret[3:5] <- quantile(x,probs=c(0.50,0.025,0.975))
    return(ret)
  }
 
  ESS = function( burnin=burnin,n_sample = n_sample){   
    #informative
    return(effectiveSize(BETA[(burnin+1):n_sample,]))        
  }
  
  
  posterior.samples <- BETA
  p = ncol(BETA)
  N = paste(c('intercept',names(diabetes)[-(which(names(diabetes)=="diabetes"))]),":") 
  colnames(posterior.samples) <-  paste(N, "beta_",0:(p-1),sep="")
  posterior.stats <- apply(posterior.samples,2,sum.stats)
  colnames(posterior.stats) <-  paste(N,"beta_",0:(p-1),sep="")
  ESS.out =round(ESS( burnin=burnin,n_sample = n_sample),0)
  ESS.out = c(ESS.out,min(ESS.out))
  names(ESS.out) =  c((paste("beta_",0:(p-1),sep="")),'minimium')
  
  ACF = matrix(NA,nrow=3,ncol=p)
  for (i in 1:p){
    ACF[,i] = acf(BETA[(burnin+1):n_sample,][,i],2,plot=FALSE)$acf
    N = paste(c('intercept',names(diabetes)[-(which(names(diabetes)=="diabetes"))]),":") 
    colnames(ACF) = paste(N, "beta_",0:(p-1),sep="")
    rownames(ACF) = paste("lag_",0:2,sep="")    
  }   
  ACF=cbind(ACF[2,],ACF[3,],100*(ACF[3,]-ACF[2,])/ACF[2,])
  colnames(ACF) = c('lag1','lag2',"% change")  
  
  return(list("posterior.stats"=round(t(posterior.stats),4),
              "posterior.samples"=posterior.samples,
              "posterior.ESS"=ESS.out,"posterior.ACF" =ACF  ))  
}

#############################FUNCTION END
n_sample = 5000
burnin = 200
MLE_logit_model = glm(diabetes ~.,family = binomial(link='logit'), data = diabetes)
b_mle =  matrix(data = summary(MLE_logit_model)$coefficients[, 1],nrow = dim(diabetes)[2]) #MLE beta 
n_beta = dim(diabetes)[2]
b_null = matrix(rep(0, n_beta),nrow = n_beta)
#############
#Sensitive Test
#############
BETA_mle_1 = gibbs_sampler(n_sample = n_sample, burnin = burnin, data = diabetes, mi=1, C=1, b0=b_mle)
BETA_mle_1000 = gibbs_sampler(n_sample = n_sample, burnin = burnin, data = diabetes,mi=1,C = 1000, b0=b_mle)
BETA_null = gibbs_sampler(n_sample = n_sample, burnin = burnin, data = diabetes, mi=1,C = 1, b0=b_null)
summary(MLE_logit_model)
confint(MLE_logit_model,level = 0.95)

BETA_mle_1$posterior.stats
BETA_mle_1000$posterior.stats
BETA_null$posterior.stats

BETA_mle_1$posterior.ESS
BETA_mle_1000$posterior.ESS
BETA_null$posterior.ESS


BETA_mle_1$posterior.ACF
BETA_mle_1000$posterior.ACF
BETA_null$posterior.ACF
ACF_Change_Compare = cbind(BETA_mle_1$posterior.ACF[,3],
      BETA_mle_1000$posterior.ACF[,3], 
      BETA_null$posterior.ACF[,3])
colnames(ACF_Change_Compare) = c("mle_1","mle_1000","null")
round(ACF_Change_Compare,1)
##############
#traceplot
#############
n_beta = dim(diabetes)[2]
#BETA1 informative prior traceplot
n_beta = dim(diabetes)[2]
par(mfrow=c(3,3))
for (i in 1:(n_beta)){
  traceplot(mcmc(BETA_mle_1$posterior.samples[,i]),
    main=paste('informative prior: beta_',(i-1)))
}
#BETA2 weak prior traceplot
par(mfrow=c(3,3))
for (i in 1:(n_beta)){
  traceplot(mcmc(BETA_mle_1000$posterior.samples[,i]),
    main=paste('weak prior: beta_',(i-1)))
}
#BETA_null prior traceplot
par(mfrow=c(3,3))
for (i in 1:(n_beta)){
  traceplot(mcmc(BETA_null$posterior.samples[,i]),
    main=paste('trivial prior: beta_',(i-1)))
}

####################################
##AFTER Parameter Selection
##Discard location, gender, BloodPressureHigh, BloodPressureLow
####################################
diabetes = read.table('Diabetes.txt', header = TRUE)[,-1]
diabetes$weight = (diabetes$hip+diabetes$weight+diabetes$waist)/3
diabetes = diabetes[-which(names(diabetes)%in%c("waist","hip","location","gender","BloodPressureHigh","BloodPressureLow"))]

n_sample = 7000
burnin = 200
M2_MLE_logit_model = glm(diabetes ~.,family = binomial(link='logit'), data = diabetes)
b_mle =  matrix(data = summary(M2_MLE_logit_model)$coefficients[, 1],nrow = dim(diabetes)[2]) #MLE beta 
n_beta = dim(diabetes)[2]
b_null = matrix(rep(0, n_beta),nrow = n_beta)
M2_BETA_null = gibbs_sampler(n_sample = n_sample, burnin = burnin, data = diabetes, mi=1,C = 1, b0=b_null)
M2_BETA_null$posterior.stats
round(exp(M2_BETA_null$posterior.stats[,1]),2)

