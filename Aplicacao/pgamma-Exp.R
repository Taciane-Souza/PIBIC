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

### distribuição de base da Gamma

G=function(t,kappa){
  return(pgamma(t, shape = kappa, scale = 1) )  
}

g=function(t,kappa){
  return(dgamma(t, shape = kappa, scale = 1) )  
}

#G=function(t,eta,kappa){
#  return(pgamma(t*eta, shape = kappa) )  
#}

#g=function(t, eta, kappa){
#  return(dgamma(t*eta, shape = kappa)*eta )  
#}


S.pGexp <- function(t,alpha,kappa,X,beta){
  
  p0 <- as.vector(exp(X %*% beta)/(1+exp(X %*% beta)))
  
  st <-  1-(1-p0)*(G(t,kappa))^alpha
  
  return(st)
}


h.pGexp=function(t,alpha,kappa,X,beta){
  p0 <- as.vector(exp(X %*% beta)/(1+exp(X %*% beta)))
  
  a=(1-p0)*alpha*g(t,kappa)*(G(t,kappa)^(alpha-1))
  
  b=1-(1-p0)*(G(t,kappa))^alpha
  return(a/b)  
}

log.lik.pGexp <- function(par,t,X,px,cens){
  kappa  <- par[1]
  alpha  <- par[2]
  nw     <- 2
  beta    <- par[(nw+1):(nw+px)]
  
  loglik <- sum( cens*log(h.pGexp(t=t,alpha=alpha,kappa=kappa,X=X,beta=beta) ) + log(S.pGexp(t=t,alpha=alpha,kappa=kappa,X=X,beta=beta)) )
  cat(-loglik,"\n")
  return(-loglik)
}

#### application ####
set.seed(123)
dados=read.table("dados_melanomaFOSP.txt",header = TRUE)
time   <- dados$tempo
cens <- dados$status
X1     <- dados$sexo 
X2     <- dados$EC_cat 
X3     <- dados$cirurgia 
X4     <- dados$radio 
X5     <- dados$quimio 
X6     <- ifelse(dados$idade>60,1,0) # ok 
X     <- model.matrix(~X1)  ## -1 --> sem beta0
px     <- ncol(X) ## number of regression parameters         

theta0=c(rep(.1,2),rep(-0.001,px)) ## initial values
length(theta0)

nrow(X)
#### estimates
estimates=ucminf(theta0, log.lik.pGexp, hessian = TRUE,t=time,X=X,px=px,cens=cens)
estimates$par

AIC = 2*estimates$value + 2*length(theta0)
AIC

BIC = 2*estimates$value + log(nrow(X))*length(theta0)
BIC

#eta     <- round(estimates$par[1],3);eta
kappa     <- round(estimates$par[1],3);kappa
alpha    <- round(estimates$par[2],3);alpha
beta     <- round(estimates$par[3:(2+px)],3);beta

# ### desvios padr?es estimados
estimates$hessian

#sd.eta<-sqrt(solve(estimates$hessian)[1,1]);round(sd.eta,3)
sd.kappa<-sqrt(solve(estimates$hessian)[1,1]);round(sd.kappa,3)
sd.alpha<-sqrt(solve(estimates$hessian)[2,2]);round(sd.alpha,3)
sd.beta0<-sqrt(solve(estimates$hessian)[3,3]);round(sd.beta0,3)
sd.beta1<-sqrt(solve(estimates$hessian)[4,4]);round(sd.beta1,3)


### intervalos de confianca de 95%
cc<-qnorm(0.975,0,1)

#ic.eta<-eta-c(1,-1)*cc*(sd.eta)
#round(ic.eta,3)


ic.kappa<-kappa-c(1,-1)*cc*(sd.kappa)
round(ic.kappa,3)


ic.alpha<-alpha-c(1,-1)*cc*(sd.alpha)
round(ic.alpha,3)

ic.beta0<-beta[1]-c(1,-1)*cc*(sd.beta0)
round(ic.beta0,3)


ic.beta1<-beta[2]-c(1,-1)*cc*(sd.beta1)
round(ic.beta1,3)


#> proporções de cura com seus ICs de 95% ####
#> # Feminino
p0F=exp(beta[1]+beta[2])/(1+exp(beta[1]+beta[2]))
p0F=round(p0F,3)
p0F
sd.p0F=deltamethod(~exp(x3+x4)/(1+exp(x3+x4)), estimates$par, cov = solve(estimates$hessian))
sd.p0F
ic.p0F<-p0F-c(1,-1)*cc*(sd.p0F)
round(ic.p0F,3)

# Masculino
p0M=exp(beta[1])/(1+exp(beta[1]))
p0M=round(p0M,3)
p0M
sd.p0M=deltamethod(~exp(x3)/(1+exp(x3)), estimates$par, cov=solve(estimates$hessian))
sd.p0M
ic.p0M<-p0M-c(1,-1)*cc*(sd.p0M)
round(ic.p0M,3)


############################
## Sobrevivências estimadas
############################
time=dados$tempo
status<-dados$status
ekm<-survfit(Surv(time,status)~X1)
time=sort(time)
st.pGexp0<-S.pGexp(t=time,alpha=alpha,kappa=kappa,X=1,beta=beta[1])
st.pGexp1<-S.pGexp(t=time,alpha=alpha,kappa=kappa,X=1,beta=beta[1]+beta[2])
plot(ekm,xlab="Tempo (anos)",conf.int=F,ylab="Sobrevivência",pch=16,lwd=2,bty="l")
lines(c(0,time),c(1,st.pGexp0),col=4,lwd=3.8)
lines(c(0,time),c(1,st.pGexp1),col=2,lwd=3.8)
#abline(h=min(st.pGexp0),lty=5,lwd=2,col=3)
#abline(h=min(st.pGexp1),lty=5,lwd=2,col=3)
legend(12,1.1, legend=c("Feminino","Masculino"),lwd = 3.8,bty = "n",col=c(2,4),cex=1.8)
legend(0,0.5, expression(hat(p)[0][F]=="0,672",IC("95%")=="[0,651; 0,693]"),pch=15,col = 2,bty = "n",cex=1.8)
legend(10,0.5, expression(hat(p)[0][M]=="0,519",IC("95%")=="[0,493; 0,545]"),pch=15,col = 4,bty = "n",cex=1.8)


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
  plot(xq2, d2s, xlab = quote("Quantis teóricos"),
       ylab = quote("Quantis empíricos"), 
       pch = 20, ylim = fy)
  par(new = T)
  plot(xq2, d21, type = "l", ylim = fy, xlab = "", ylab = "", lwd=1.2)
  par(new = T)
  plot(xq2, d2med, type = "l", ylim = fy, xlab = "", ylab = "", lwd=2)
  par(new = T)
  plot(xq2, d22, type = "l", ylim = fy, xlab = "", ylab = "", lwd=1.2)
}


## residuals plot 
sobre=S.pGexp(t=dados$tempo,alpha=alpha,kappa=kappa,X=X,beta=beta)
qr <- qnorm(dados$status*(1 - sobre) + (1-dados$status)*runif(length(time),1-sobre))
envelopeDS(qr) 




## outro tipo de resíduo
# ### residuos de cox-snell
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

