rm(list = ls()) # limpa os dados da memória do R

##############################################################
## Packages
##############################################################

library(numDeriv)
library(labstatR)
library(MASS)
library(robustbase)
library(maxLik) 
library(moments)
library(zipfR)
library(msm)
#library(optimr)
library(ucminf)
library(optimx)
library(survival)

### distribuição de base da Weibull
G=function(t,eta,kappa){
  return(1-exp(-kappa*t^(eta)) )  
}


g=function(t,eta,kappa){
  return(kappa*eta*t^(eta-1)*exp(-kappa*t^eta) )  
}

st.mix <- function(t,eta,kappa,X,beta){
  
  p <- as.vector(exp(X %*% beta)/(1+exp(X %*% beta)))
  
  st <-  p+(1-p)*(1-G(t,eta,kappa))
  
  return(st)
}


ht.mix=function(t,eta,kappa,X,beta){
  p <- as.vector(exp(X %*% beta)/(1+exp(X %*% beta)))
  
  a=(1-p)*g(t,eta,kappa)
  
  b=st.mix(t,eta,kappa,X,beta)
  return(a/b)  
}


log.lik.mix <- function(par,t,X,px,cens){
  eta  <- par[1]
  kappa  <- par[2]
  nw     <- 2
  beta    <- par[(nw+1):(nw+px)]
  
  loglik <- sum( cens*log(ht.mix(t=t,eta=eta,kappa=kappa,X=X,beta=beta) ) + log(st.mix(t=t,eta=eta,kappa=kappa,X=X,beta=beta)) )
  cat(-loglik,"\n")
  return(-loglik)
}


#### application ####
setwd("G:/My Drive/PIBIC-Taciane/Codigos/application")
dados=read.table("dados_melanomaFOSP.txt",header = TRUE)
time   <- dados$tempos
cens <- dados$status
X1     <- dados$sexo
X2     <- dados$EC_cat 
X3     <- dados$cirurgia 
X4     <- dados$radio 
X5     <- dados$quimio 
X6     <- ifelse(dados$idade>60,1,0) # ok 
X     <- model.matrix(~X1)  ## -1 --> sem beta0
px     <- ncol(X) ## number of regression parameters         

p0=c(rep(1,2),rep(-0.001,px)) ## initial values
length(p0)


#### estimates
estimates=ucminf(p0, log.lik.mix, hessian = TRUE,t=time,X=X,px=px,cens=cens)
estimates$par

AIC = 2*estimates$value + 2*length(p0)
AIC

eta     <- round(estimates$par[1],3);eta
kappa     <- round(estimates$par[2],3);kappa
beta     <- round(estimates$par[3:(2+px)],3);beta

# ### desvios padr?es estimados

sd.eta<-sqrt(solve(estimates$hessian)[1,1]);round(sd.eta,3)

sd.kappa<-sqrt(solve(estimates$hessian)[2,2]);round(sd.kappa,3)
sd.beta0<-sqrt(solve(estimates$hessian)[3,3]);round(sd.beta0,3)
sd.beta1<-sqrt(solve(estimates$hessian)[4,4]);round(sd.beta1,3)
sd.beta2<-sqrt(solve(estimates$hessian)[5,5]);round(sd.beta2,3)

### intervalos de confianca de 95%
cc<-qnorm(0.975,0,1)

ic.eta<-eta-c(1,-1)*cc*(sd.eta)
round(ic.eta,3)


ic.kappa<-kappa-c(1,-1)*cc*(sd.kappa)
round(ic.kappa,3)


ic.beta0<-beta[1]-c(1,-1)*cc*(sd.beta0)
round(ic.beta0,3)


ic.beta1<-beta[2]-c(1,-1)*cc*(sd.beta1)
round(ic.beta1,3)

ic.beta2<-beta[3]-c(1,-1)*cc*(sd.beta2)
round(ic.beta2,3)

# par(mfrow=c(1,2))
# ############################
# ## Sobrevivências estimadas
# ############################
# 
# time=dados$tempo_anos
# status<-dados$status
# ekm<-survfit(Surv(time,status)~X6)
# time=sort(time)
# st.pGexp0<-S.pGexp(t=time,alpha=alpha,eta=eta,kappa=kappa,X=1,beta=beta[1])
# st.pGexp1<-S.pGexp(t=time,alpha=alpha,eta=eta,kappa=kappa,X=1,beta=beta[1]+beta[2])
# plot(ekm,xlab="Tempo (anos)",conf.int=F,ylab="Survival",pch=16,lwd=2,bty="l")
# lines(c(0,time),c(1,st.pGexp0),col=4,lwd=2,lty=3)
# lines(c(0,time),c(1,st.pGexp1),col=2,lwd=2,lty=2)
# abline(h=min(st.pGexp0),lty=5,lwd=2,col=3)
# abline(h=min(st.pGexp1),lty=5,lwd=2,col=3)
# legend(10,1.0, legend=c("<=60 anos",">60 anos"),lty = c(2,3),lwd = 3,bty = "n",col=c(2,4))


#### Residual analysis
### residual plot + envelope
## argument x=residuals
envelopeDS <- function(x){
  U	         <- x
  n	         <- length(x)
  d2s 	     <- sort(U)
  xq2 	     <- qnorm(ppoints(n))
  Xsim 	     <- matrix(0, 100, n)
  for(i in 1:100){
    u2       <- rnorm(n)
    Xsim[i,] <- u2
  }
  Xsim2      <- apply(Xsim, 1, sort)
  d21        <- matrix(0, n, 1)
  d22        <- matrix(0, n, 1)
  for(i in 1:n){
    d21[i]  <- quantile(Xsim2[i,], 0.025)
    d22[i]  <- quantile(Xsim2[i,], 0.975)
  }
  d2med      <- apply(Xsim2, 1, mean)
  fy         <- range(d2s, d21, d22)
  plot(xq2, d2s, xlab = quote("Theoretical Quantiles"),
       ylab = quote("Quantile Residuals"), 
       pch = 20, ylim = fy)
  par(new = T)
  plot(xq2, d21, type = "l", ylim = fy, xlab = "", ylab = "", lwd=1.2)
  par(new = T)
  plot(xq2, d2med, type = "l", ylim = fy, xlab = "", ylab = "", lwd=2)
  par(new = T)
  plot(xq2, d22, type = "l", ylim = fy, xlab = "", ylab = "", lwd=1.2)
}


## residuals plot 

sobre=st.mix(t=dados$tempo_anos,eta=eta,kappa=kappa,X=X,beta=beta)
qr <- qnorm(dados$status*(1 - sobre) + (1-dados$status)*runif(length(time),1-sobre))
envelopeDS(qr) 


# # outro tipo de resíduo
# ### residuos de cox-snell
# windows()
# library(ReIns)
# ei    = -log(S.pGexp(t=dados$tempo_anos,alpha=alpha,eta=eta,kappa=kappa,X=X,beta=beta))
# censored=dados$status==0
# censored
# 
# a=cExpQQ(ei, censored, plot = FALSE, main = "Exponential QQ-plot")
# 
# par(mfrow=c(1,2))
# plot(a$eqq.the,a$eqq.emp,xlab="Theoretical quantiles",ylab="Empirical quantiles",pch=16,col=1,)
# abline(c(0,1),c(1,0),col=1,lwd=2)
# 


p00=as.vector(exp(X %*% beta)/(1+exp(X %*% beta)))
prop.cura=mean(p00)
prop.cura








