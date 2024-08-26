rm(list=ls()) 
if(!require(survival)){install.packages("survival")}
if(!require(tcltk2)){install.packages("tcltk2")}
if(!require(numDeriv)){install.packages("numDeriv")}
if(!require(xlsx)){install.packages("xlsx")}
if(!require(msm)){install.packages("msm")}
 
  t0=Sys.time()
  name="pWeibull-Exp"
  # Os cenários foram criados com base nos seguintes parâmetros:
  # eta, kappa, alpha, β0, β1
  
  st1=c(1, .5, 1,   1, 1)    # 1º cenário
  st2=c(1,  1, 1,   1, 1)    # 2º cenário
  st3=c(1,1.5, 1.5, 1, 1)    # 3º cenário
  setups=rbind(st1,st2,st3)

  
  S=2        # número de replicações
  ni=100     # tamanho inicial da amostra
  nf=2500    # tamanho final da amostra
  pass=200   # salto entre cada amostra
  samples=seq(ni,nf,pass)

  den=function(t,e,k){ dweibull(t,k,(1/e)^(1/k) ) } # função densidade da weibull
  
  cum=function(t,e,k){ pweibull(t,k,(1/e)^(1/k)) } # função acumulada da weibull

  # a = alpha; e = eta; k = kappa.
  # densidade da pweibull
  propden=function(a,e,k,b0,b1,x,t){
  p=exp(b0+b1*x)/(1+exp(b0+b1*x))
  f0=den(t,e,k)
  c0=cum(t,e,k)
  cc=(1-p)^(1/a)*a*f0*(((1-p)^(1/a))*c0)^(a-1)
  return(cc)
  }
  
  # acumulada da pweibull
  propcum=function(a,e,k,b0,b1,x,t){
  p=exp(b0+b1*x)/(1+exp(b0+b1*x))
  f0=den(t,e,k)
  c0=cum(t,e,k)
  dd=(((1-p)^(1/a))*c0)^(a)
  return(dd)
  }

  # função de fração de cura
  cure<-function(a,e,k,b0,b1,x){
    p=exp(b0+b1*x)/(1+exp(b0+b1*x))
    ee=(1-p)^(1/a)
    return(ee)
    }
  
  #  Auxilia na geração dos valores
  f<-function(a,e,k,b0,b1,x,unif,t){propcum(a,e,k,b0,b1,x,t)-unif}

  # função para gerar um determinante
  # - para cada parâmetro calcula-se uma função que a partir dela serão geradas 
  # as matrizes jacobianas para cada um dos parâmetros
  #est[1] = estimativa para eta; est[2] = estimativa para kappa; est[4] = b0; est[5] = b1
  md=function(est,covmat,x){
  myfun.e=function(e){1-pweibull(1/(1+exp(-est[4]-est[5]*x)) ,shape=e,scale=est[2])}
  myfun.k=function(k){1-pweibull(1/(1+exp(-est[4]-est[5]*x)),shape=est[1],scale=k)}
  myfun.a=function(a){1-pweibull(1/(1+exp(-est[4]-est[5]*x)),shape=est[1],scale=est[2])}
  myfun.b0=function(b0){1-pweibull(1/(1+exp(-b0-est[5]*x)),shape=est[1],scale=est[2])}
  myfun.b1=function(b1){1-pweibull(1/(1+exp(-est[4]-b1*x)),shape=est[1],scale=est[2])}
  j1=jacobian(myfun.e,est[1])
  j2=jacobian(myfun.k,est[2])
  j3=jacobian(myfun.a,est[3])
  j4=jacobian(myfun.b0,est[4])
  j5=jacobian(myfun.b1,est[5])
  j=c(j1,j2,j3,j4,j5) # vetor de matrizes
  v=j%*%covmat%*%j
  se=sqrt(v)
  return(se)}

  # função log-verossimilhança
  logvero<-function(par){
  e=par[1]
  k=par[2]
  a=par[3]
  b0=par[4]
  b1=par[5]
  f=propden(a,e,k,b0,b1,x1,t1)  #  ???
  s=1-propcum(a,e,k,b0,b1,x1,t1)
  logl=sum((c1*log(f))+((1-c1)*log(s)))
  logl=ifelse(logl=="NaN",-10^10,logl)
  return(-logl)}
  
  # matriz de uns com 3 linhas e 2 colunas
  p0s=matrix(1,nrow(setups),2)
  # setups é a matriz contendo os três cenários
  
  #> este *for* termina na linha 339
  for(setup in 1:nrow(setups)){

  e0=setups[setup,1]
  k0=setups[setup,2]
  a0=setups[setup,3]
  b00=setups[setup,4]
  b10=setups[setup,5]

  p00=cure(a0,e0,k0,b00,b10,0)
  p10=cure(a0,e0,k0,b00,b10,1)

  p0s[setup,1]=p00
  p0s[setup,2]=p10

  st=c(a0,e0,k0,b00,b10,p00,p10)




    cover.a=numeric(); cover.e=numeric();  cover.k=numeric(); cover.b0=numeric(); cover.b1=numeric(); cover.p0=numeric(); cover.p1=numeric()
      mse.a=numeric();   mse.e=numeric();    mse.k=numeric();   mse.b0=numeric();   mse.b1=numeric();   mse.p0=numeric();   mse.p1=numeric()
     bias.a=numeric();  bias.e=numeric();   bias.k=numeric();  bias.b0=numeric();  bias.b1=numeric();  bias.p0=numeric();  bias.p1=numeric()
       ic.a=numeric();    ic.e=numeric();     ic.k=numeric();    ic.b0=numeric();    ic.b1=numeric();    ic.p0=numeric();    ic.p1=numeric()




  error=numeric()
  censor=numeric()
  ind=0

  for(n in samples){
  ind=ind+1
  pb <- tkProgressBar(title = "Simulation Progress", min = 0, max = S, width = 1000)

  a=numeric()
  e=numeric()
  k=numeric()
  be0=numeric()
  be1=numeric()
  p0=numeric()
  p1=numeric()

  da=numeric()
  de=numeric()
  dk=numeric()
  db0=numeric()
  db1=numeric()
  dp0=numeric()
  dp1=numeric()

  cont_a=0
  cont_e=0
  cont_k=0
  cont_b0=0
  cont_b1=0
  cont_p0=0
  cont_p1=0

  avg=0
  eta=0
  errors=0
  counter=0
  m=1
  times=numeric()
  cens=numeric()

  while(m <= S){
  counter=counter+1

  setTkProgressBar(pb, m-1, label=paste("Setup:",setup,"/",nrow(setups),"         ","Sample size:",n,"(",ind,"/",length(samples),")",
  "         ","Simulations done:",m-1,"/",S,"         ","Errors:",errors,"         ","ATS:", avg,"secs","         ","ETA:",eta,"mins"))

  ti=Sys.time()

  #Data generation
    t1=numeric()
    c1=numeric()
    x1=rbinom(n,1,0.5)
    tp=numeric()
    tk=numeric()

    for(i in 1:n){
      try(
        if(x1[i]==0){ if(rbinom(1,1,(1-p00)^(1/a0))==0){tp[i]=Inf}
          else {teste=uniroot(f,c(0,10^6),tol=10^(-6), a=a0,e=e0,k=k0,b0=b00,b1=b10,x=x1[i],unif=runif(1,0,(1-p00)^(1/a0))); tp[i]=teste$root}}
        else { if(rbinom(1,1,(1-p10)^(1/a0))==0){tp[i]=Inf}
          else {teste=uniroot(f,c(0,10^6),tol=10^(-6), a=a0,e=e0,k=k0,b0=b00,b1=b10,x=x1[i],unif=runif(1,0,(1-p10)^(1/a0))); tp[i]=teste$root}}
      , silent = TRUE)
      }
    

    for(i in 1:length(tp)){ if(tp[i]==Inf | is.na(tp[i]==T)){ tk[i]=0} else tk[i]=tp[i] }
    cen=runif(n,0,max(tk))
    for(j in 1:n){t1[j]=min(tp[j],cen[j]); if(tp[j]<=cen[j] | is.na(tp[j]==T)) c1[j]=1 else c1[j]=0}
    cens[m]=1-sum(c1)/n

  ekm = survfit(Surv(t1,c1)~x1, se.fit = FALSE)
  plot(ekm)

  #MLE Estimator
    mod<- try(optim(c(a0,e0,k0,b00,b10), logvero, method="BFGS", hessian=T), silent=T);
    if( is(mod,"try-error")==T ){errors = errors+1; next}
    if( mod$convergence!=0 ){errors = errors+1; next}
    summary(mod)
    cov=try(solve(mod$hessian), silent=T);
    if( is(cov,"try-error")==T ){errors = errors+1; next}
  #Keeping the values of interest
    if( mod$par[1]<=0 | mod$par[2]<=0 | mod$par[3]<=0 |
        is.na(cov[1,1])==T | is.na(cov[2,2])==T | is.na(cov[3,3])==T | is.na(cov[4,4])==T | is.na(cov[5,5])==T |
        cov[1,1]<0 | cov[2,2]<0 | cov[3,3]<0 | cov[4,4]<0 | cov[5,5]<0 ){errors=errors+1; next}
    else {
    e1=mod$par[1]
    k1=mod$par[2]
    a1=mod$par[3]
    b01=mod$par[4]
    b11=mod$par[5]

    a[m]=a1
    e[m]=e1
    k[m]=k1
    be0[m]=b01
    be1[m]=b11
    p0[m]=cure(a[m],e[m],k[m],be0[m],be1[m],0)
    p1[m]=cure(a[m],e[m],k[m],be0[m],be1[m],1)

    de[m]=sqrt(cov[1,1])
    dk[m]=sqrt(cov[2,2])
    da[m]=sqrt(cov[3,3])
    db0[m]=sqrt(cov[4,4])
    db1[m]=sqrt(cov[5,5])
    dp0[m]=(deltamethod (~(1-(exp(x4)/(1+exp(x4)) ))^(1/x1), mod$par, solve(mod$hessian)))^2; if( dp0[m]=="NaN" ){errors = errors+1; next}
    dp1[m]=(deltamethod (~(1-(exp(x4+x5)/(1+exp(x4+x5)) ))^(1/x1), mod$par, solve(mod$hessian)))^2; if( dp1[m]=="NaN" ){errors = errors+1; next}

    int_a=a[m]+c(-1,1)*qnorm(0.975)*da[m]
    int_e=e[m]+c(-1,1)*qnorm(0.975)*de[m]
    int_k=k[m]+c(-1,1)*qnorm(0.975)*dk[m]
    int_b0=be0[m]+c(-1,1)*qnorm(0.975)*db0[m]
    int_b1=be1[m]+c(-1,1)*qnorm(0.975)*db1[m]
    int_p0=p0[m]+c(-1,1)*qnorm(0.975)*dp0[m]
    int_p1=p1[m]+c(-1,1)*qnorm(0.975)*dp1[m]

    if(a0>=int_a[1] & a0<=int_a[2]) cont_a=cont_a+1
    if(e0>=int_e[1] & e0<=int_e[2]) cont_e=cont_e+1
    if(k0>=int_k[1] & k0<=int_k[2]) cont_k=cont_k+1
    if(b00>=int_b0[1] & b00<=int_b0[2]) cont_b0=cont_b0+1
    if(b10>=int_b1[1] & b10<=int_b1[2]) cont_b1=cont_b1+1
    if(p00>=int_p0[1] & p00<=int_p0[2]) cont_p0=cont_p0+1
    if(p10>=int_p1[1] & p10<=int_p1[2]) cont_p1=cont_p1+1


    tf=Sys.time()
    times[m]=tf-ti
    avg=round(mean(times),2)
    eta=round((avg*(S-m+(errors/m)*(S-m)))/60,2)
    m=m+1
    }

    }
    close(pb)
    error[ind]=errors
    censor[ind]=mean(cens)


    cover.a[ind]=cont_a/S
    cover.e[ind]=cont_e/S
    cover.k[ind]=cont_k/S
    cover.b0[ind]=cont_b0/S
    cover.b1[ind]=cont_b1/S
    cover.p0[ind]=cont_p0/S
    cover.p1[ind]=cont_p1/S
    mse.a[ind]=(mean(a)-a0)^2+mean(da)^2
    mse.e[ind]=(mean(e)-e0)^2+mean(de)^2
    mse.k[ind]=(mean(k)-k0)^2+mean(dk)^2
    mse.b0[ind]=(mean(be0)-b00)^2+mean(db0)^2
    mse.b1[ind]=(mean(be1)-b10)^2+mean(db1)^2
    mse.p0[ind]=(mean(p0)-p00)^2+mean(dp0)^2
    mse.p1[ind]=(mean(p1)-p10)^2+mean(dp1)^2
    bias.a[ind]=mean(a)-a0
    bias.e[ind]=mean(e)-e0
    bias.k[ind]=mean(k)-k0
    bias.b0[ind]=mean(be0)-b00
    bias.b1[ind]=mean(be1)-b10
    bias.p0[ind]=mean(p0)-p00
    bias.p1[ind]=mean(p1)-p10
    ic.a[ind]=2*qnorm(0.975)*mean(da)
    ic.e[ind]=2*qnorm(0.975)*mean(de)
    ic.k[ind]=2*qnorm(0.975)*mean(dk)
    ic.b0[ind]=2*qnorm(0.975)*mean(db0)
    ic.b1[ind]=2*qnorm(0.975)*mean(db1)
    ic.p0[ind]=2*qnorm(0.975)*mean(dp0)
    ic.p1[ind]=2*qnorm(0.975)*mean(dp1)

  }

  a=bias.a+a0
  e=bias.e+e0
  k=bias.k+k0
  b0=bias.b0+b00
  b1=bias.b1+b10
  p0=bias.p0+p00
  p1=bias.p1+p10

  da=ic.a/(2*qnorm(0.975))
  de=ic.e/(2*qnorm(0.975))
  dk=ic.k/(2*qnorm(0.975))
  db0=ic.b0/(2*qnorm(0.975))
  db1=ic.b1/(2*qnorm(0.975))
  dp0=ic.p0/(2*qnorm(0.975))
  dp1=ic.p1/(2*qnorm(0.975))


  tab=round(cbind(error,censor,samples,a,e,k,b0,b1,p0,p1,da,de,dk,db0,db1,dp0,dp1),4)
  aux1=rep("",ncol(tab))
  aux2=c("Initial Values:","","",a0,st[2],k0,b00,b10,round(p00,4),round(p10,4),"","","","","Simulations:","",S)
  tabela=rbind(tab,aux1,aux2)




  #Saving the good stuff
  write.xlsx(tabela, paste("simulation.",name,".xlsx", sep=""), paste("Setup",setup, sep=""), row.names=F, append=T )



assign(paste("cover.a.", setup, sep = ""), cover.a)
assign(paste("cover.e.", setup, sep = ""), cover.e)
assign(paste("cover.k.", setup, sep = ""), cover.k)
assign(paste("cover.b0.", setup, sep = ""), cover.b0)
assign(paste("cover.b1.", setup, sep = ""), cover.b1)
assign(paste("cover.p0.", setup, sep = ""), cover.p0)
assign(paste("cover.p1.", setup, sep = ""), cover.p1)

assign(paste("mse.a.", setup, sep = ""), mse.a)
assign(paste("mse.e.", setup, sep = ""), mse.e)
assign(paste("mse.k.", setup, sep = ""), mse.k)
assign(paste("mse.b0.", setup, sep = ""), mse.b0)
assign(paste("mse.b1.", setup, sep = ""), mse.b1)
assign(paste("mse.p0.", setup, sep = ""), mse.p0)
assign(paste("mse.p1.", setup, sep = ""), mse.p1)

assign(paste("bias.a.", setup, sep = ""), bias.a)
assign(paste("bias.e.", setup, sep = ""), bias.e)
assign(paste("bias.k.", setup, sep = ""), bias.k)
assign(paste("bias.b0.", setup, sep = ""), bias.b0)
assign(paste("bias.b1.", setup, sep = ""), bias.b1)
assign(paste("bias.p0.", setup, sep = ""), bias.p0)
assign(paste("bias.p1.", setup, sep = ""), bias.p1)

assign(paste("ic.a.", setup, sep = ""), ic.a)
assign(paste("ic.e.", setup, sep = ""), ic.e)
assign(paste("ic.k.", setup, sep = ""), ic.k)
assign(paste("ic.b0.", setup, sep = ""), ic.b0)
assign(paste("ic.b1.", setup, sep = ""), ic.b1)
assign(paste("ic.p0.", setup, sep = ""), ic.p0)
assign(paste("ic.p1.", setup, sep = ""), ic.p1)

}

#> Gráficos ####
col1=2
col2=3
col3=4
pch1=15
pch2=16
pch3=17

      pdf(paste("simulation.",name,".pdf", sep=""),width=10, height=14)

      par(mfrow=c(7,4), mar=c(1.2,4,2,0.5), oma=c(1,1,1,0))

      plot(samples,mse.a.1,type="o",pch=pch1, col=col1, ylab="Parameter alpha", ylim=c(0,max(mse.a.1,mse.a.2,mse.a.3)))
      lines(samples,mse.a.2, type="o", pch=pch2, col=col2 )
      lines(samples,mse.a.3, type="o", pch=pch3, col=col3 )
      abline(h=0,lty=2)
      title("Mean Square Errors", cex.main = 1,   font.main= 1)
      plot(samples,bias.a.1,type="o",pch=pch1, col=col1,ylab="", ylim=c(min(bias.a.1,bias.a.2,bias.a.3),max(bias.a.1,bias.a.2,bias.a.3)))
      lines(samples,bias.a.2, type="o", pch=pch2, col=col2 )
      lines(samples,bias.a.3, type="o", pch=pch3, col=col3 )
      abline(h=0,lty=2)
      title("Biases", cex.main = 1,   font.main= 1)
      plot(samples,cover.a.1,type="o",pch=pch1, col=col1,ylab="", xlab="", ylim=c(0.9,1.0))
      lines(samples,cover.a.2, type="o", pch=pch2, col=col2 )
      lines(samples,cover.a.3, type="o", pch=pch3, col=col3 )
      abline(h=0.95,lty=2)
      abline(h=0.936,lty=3)
      abline(h=0.964,lty=3)
      title("Coverage Probabilities", cex.main = 1,   font.main= 1)
      plot(samples,ic.a.1,type="o",pch=pch1, col=col1,ylab="", ylim=c(0,max(ic.a.1,ic.a.2,ic.a.3)))
      lines(samples,ic.a.2, type="o", pch=pch2, col=col2 )
      lines(samples,ic.a.3, type="o", pch=pch3, col=col3 )
      abline(h=0,lty=2)
      title("Coverage Lengths", cex.main = 1,   font.main= 1)

      plot(samples,mse.e.1,type="o",pch=pch1, col=col1, ylab="Parameter eta", ylim=c(0,max(mse.e.1,mse.e.2,mse.e.3)))
      lines(samples,mse.e.2, type="o", pch=pch2, col=col2 )
      lines(samples,mse.e.3, type="o", pch=pch3, col=col3 )
      abline(h=0,lty=2)
      plot(samples,bias.e.1,type="o",pch=pch1, col=col1, ylab="", ylim=c(min(bias.e.1,bias.e.2,bias.e.3),max(bias.e.1,bias.e.2,bias.e.3)))
      lines(samples,bias.e.2, type="o", pch=pch2, col=col2 )
      lines(samples,bias.e.3, type="o", pch=pch3, col=col3 )
      abline(h=0,lty=2)
      plot(samples,cover.e.1,type="o",pch=pch1, col=col1, ylim=c(0.9,1.0), ylab="")
      lines(samples,cover.e.2, type="o", pch=pch2, col=col2 )
      lines(samples,cover.e.3, type="o", pch=pch3, col=col3 )
      abline(h=0.95,lty=2)
      abline(h=0.936,lty=3)
      abline(h=0.964,lty=3)
      plot(samples,ic.e.1,type="o",pch=pch1, col=col1,ylab="", ylim=c(0,max(ic.e.1,ic.e.2,ic.e.3)))
      lines(samples,ic.e.2, type="o", pch=pch2, col=col2 )
      lines(samples,ic.e.3, type="o", pch=pch3, col=col3 )
      abline(h=0,lty=2)

      plot(samples,mse.k.1,type="o",pch=pch1, col=col1, ylab=expression('Parameter kappa'), ylim=c(0,max(mse.k.1,mse.k.2,mse.k.3)))
      lines(samples,mse.k.2, type="o", pch=pch2, col=col2 )
      lines(samples,mse.k.3, type="o", pch=pch3, col=col3 )
      abline(h=0,lty=2)
      title("", cex.main = 1,   font.main= 1)
      plot(samples,bias.k.1,type="o",pch=pch1, col=col1, ylab="", ylim=c(min(bias.k.1,bias.k.2,bias.k.3),max(bias.k.1,bias.k.2,bias.k.3)))
      lines(samples,bias.k.2, type="o", pch=pch2, col=col2 )
      lines(samples,bias.k.3, type="o", pch=pch3, col=col3 )
      abline(h=0,lty=2)
      plot(samples,cover.k.1,type="o",pch=pch1, col=col1, ylim=c(0.9,1.0), ylab="")
      lines(samples,cover.k.2, type="o", pch=pch2, col=col2 )
      lines(samples,cover.k.3, type="o", pch=pch3, col=col3 )
      abline(h=0.95,lty=2)
      abline(h=0.936,lty=3)
      abline(h=0.964,lty=3)
      plot(samples,ic.k.1,type="o",pch=pch1, col=col1,ylab="", ylim=c(0,max(ic.k.1,ic.k.2,ic.k.3)))
      lines(samples,ic.k.2, type="o", pch=pch2, col=col2 )
      lines(samples,ic.k.3, type="o", pch=pch3, col=col3 )
      abline(h=0,lty=2)

      plot(samples,mse.b0.1,type="o",pch=pch1, col=col1, ylab=expression('Parameter '~beta[0]), ylim=c(0,max(mse.b0.1,mse.b0.2,mse.b0.3)))
      lines(samples,mse.b0.2, type="o", pch=pch2, col=col2 )
      lines(samples,mse.b0.3, type="o", pch=pch3, col=col3 )
      abline(h=0,lty=2)
      title("", cex.main = 1,   font.main= 1)
      plot(samples,bias.b0.1,type="o",pch=pch1, col=col1, ylab="", ylim=c(min(bias.b0.1,bias.b0.2,bias.b0.3),max(bias.b0.1,bias.b0.2,bias.b0.3)))
      lines(samples,bias.b0.2, type="o", pch=pch2, col=col2 )
      lines(samples,bias.b0.3, type="o", pch=pch3, col=col3 )
      abline(h=0,lty=2)
      plot(samples,cover.b0.1,type="o",pch=pch1, col=col1, ylim=c(0.9,1.0), ylab="")
      lines(samples,cover.b0.2, type="o", pch=pch2, col=col2 )
      lines(samples,cover.b0.3, type="o", pch=pch3, col=col3 )
      abline(h=0.95,lty=2)
      abline(h=0.936,lty=3)
      abline(h=0.964,lty=3)
      plot(samples,ic.b0.1,type="o",pch=pch1, col=col1,ylab="", ylim=c(0,max(ic.b0.1,ic.b0.2,ic.b0.3)))
      lines(samples,ic.b0.2, type="o", pch=pch2, col=col2 )
      lines(samples,ic.b0.3, type="o", pch=pch3, col=col3 )
      abline(h=0,lty=2)

      plot(samples,mse.b1.1,type="o",pch=pch1, col=col1, ylab=expression('Parameter '~beta[1]), ylim=c(0,max(mse.b1.1,mse.b1.2,mse.b1.3)))
      lines(samples,mse.b1.2, type="o", pch=pch2, col=col2 )
      lines(samples,mse.b1.3, type="o", pch=pch3, col=col3 )
      abline(h=0,lty=2)
      title("", cex.main = 1,   font.main= 1)
      plot(samples,bias.b1.1,type="o",pch=pch1, col=col1, ylab="", ylim=c(min(bias.b1.1,bias.b1.2,bias.b1.3),max(bias.b1.1,bias.b1.2,bias.b1.3)))
      lines(samples,bias.b1.2, type="o", pch=pch2, col=col2 )
      lines(samples,bias.b1.3, type="o", pch=pch3, col=col3 )
      abline(h=0,lty=2)
      plot(samples,cover.b1.1,type="o",pch=pch1, col=col1, ylim=c(0.9,1.0), ylab="")
      lines(samples,cover.b1.2, type="o", pch=pch2, col=col2 )
      lines(samples,cover.b1.3, type="o", pch=pch3, col=col3 )
      abline(h=0.95,lty=2)
      abline(h=0.936,lty=3)
      abline(h=0.964,lty=3)
      plot(samples,ic.b1.1,type="o",pch=pch1, col=col1,ylab="", ylim=c(0,max(ic.b1.1,ic.b1.2,ic.b1.3)))
      lines(samples,ic.b1.2, type="o", pch=pch2, col=col2 )
      lines(samples,ic.b1.3, type="o", pch=pch3, col=col3 )
      abline(h=0,lty=2)

      plot(samples,mse.p0.1,type="o",pch=pch1, col=col1, ylab=expression('Parameter '~p[0]), ylim=c(0,max(mse.p0.1,mse.p0.2,mse.p0.3)))
      lines(samples,mse.p0.2, type="o", pch=pch2, col=col2 )
      lines(samples,mse.p0.3, type="o", pch=pch3, col=col3 )
      abline(h=0,lty=2)
      title("", cex.main = 1,   font.main= 1)
      plot(samples,bias.p0.1,type="o",pch=pch1, col=col1, ylab="", ylim=c(min(bias.p0.1,bias.p0.2,bias.p0.3),max(bias.p0.1,bias.p0.2,bias.p0.3)))
      lines(samples,bias.p0.2, type="o", pch=pch2, col=col2 )
      lines(samples,bias.p0.3, type="o", pch=pch3, col=col3 )
      abline(h=0,lty=2)
      plot(samples,cover.p0.1,type="o",pch=pch1, col=col1, ylim=c(0.9,1.0), ylab="")
      lines(samples,cover.p0.2, type="o", pch=pch2, col=col2 )
      lines(samples,cover.p0.3, type="o", pch=pch3, col=col3 )
      abline(h=0.95,lty=2)
      abline(h=0.936,lty=3)
      abline(h=0.964,lty=3)
      plot(samples,ic.p0.1,type="o",pch=pch1, col=col1,ylab="", ylim=c(0,max(ic.p0.1,ic.p0.2,ic.p0.3)))
      lines(samples,ic.p0.2, type="o", pch=pch2, col=col2 )
      lines(samples,ic.p0.3, type="o", pch=pch3, col=col3 )
      abline(h=0,lty=2)

      plot(samples,mse.p1.1,type="o",pch=pch1, col=col1, ylab=expression('Parameter '~p[1]), ylim=c(0,max(mse.p1.1,mse.p1.2,mse.p1.3)))
      lines(samples,mse.p1.2, type="o", pch=pch2, col=col2 )
      lines(samples,mse.p1.3, type="o", pch=pch3, col=col3 )
      abline(h=0,lty=2)
      title("", cex.main = 1,   font.main= 1)
      plot(samples,bias.p1.1,type="o",pch=pch1, col=col1, ylab="", ylim=c(min(bias.p1.1,bias.p1.2,bias.p1.3),max(bias.p1.1,bias.p1.2,bias.p1.3)))
      lines(samples,bias.p1.2, type="o", pch=pch2, col=col2 )
      lines(samples,bias.p1.3, type="o", pch=pch3, col=col3 )
      abline(h=0,lty=2)
      plot(samples,cover.p1.1,type="o",pch=pch1, col=col1, ylim=c(0.9,1.0), ylab="")
      lines(samples,cover.p1.2, type="o", pch=pch2, col=col2 )
      lines(samples,cover.p1.3, type="o", pch=pch3, col=col3 )
      abline(h=0.95,lty=2)
      abline(h=0.936,lty=3)
      abline(h=0.964,lty=3)
      plot(samples,ic.p1.1,type="o",pch=pch1, col=col1,ylab="", ylim=c(0,max(ic.p1.1,ic.p1.2,ic.p1.3)))
      lines(samples,ic.p1.2, type="o", pch=pch2, col=col2 )
      lines(samples,ic.p1.3, type="o", pch=pch3, col=col3 )
      abline(h=0,lty=2)

      dev.off()



save.image(paste("Simulação",name,".Rdata", sep=""))
setups=cbind(setups,round(p0s,4));
colnames(setups)=c("a","e","k","beta0","beta1","p0","p1")
setups

l0=Sys.time()-t0
l1=(Sys.time()-t0)*1000/S
list("Simulation time"=l0,"Estimated time for 1000 simulations"=as.numeric(l1,unit="hours"))




