rm(list=ls(all=TRUE))
require(coda)
require(maxLik)
require(VGAM)
require(MCMCpack)
############ Declarando as variaveis e alguns parametros ##############
old  =  options(digits=20)
pphi <- 1   
pvarp <- 2  
B <- 100        # Qtd de estimativas que ser?o calculados
nini<-5
nmax<-15
jumpn<-5
lim<-50
NN<-(nmax-nini)/jumpn +1
cc<-qnorm(0.975,0,1)
##Parametros Bayesiano###
R<-5500         #Numero de iteracoes
burnin<-500
jump<-5
set.seed(2023)
ES<-4
######### tamanho das matrizes que armazenarao estimativas######
mrephi <- matrix(nrow=NN,ncol=ES)
msephi <- matrix(nrow=NN,ncol=ES)
cobphi <- matrix(nrow=NN,ncol=ES)
mrevarp <- matrix(nrow=NN,ncol=ES)
msevarp <- matrix(nrow=NN,ncol=ES)
cobvarp <- matrix(nrow=NN,ncol=ES)
gematriz1 <- matrix(nrow=NN,ncol=B)
gematriz2 <- matrix(nrow=NN,ncol=B)
gematriz3 <- matrix(nrow=NN,ncol=B)
gematriz4 <- matrix(nrow=NN,ncol=B)
gematriz5 <- matrix(nrow=NN,ncol=B)
gematriz6 <- matrix(nrow=NN,ncol=B)

posterior1 <- function (v) {
  p<- (-n*media/v)-(n+1)*log(v)-n*log(sum(exp(-x/v)))
  return(p) }

posterior2 <- function (v,phi) {
  p<--(n+2)*log(v)-(n*(media-phi)/v)-sum(exp(-(x-phi)/v))
  return(p)
}

posterior3 <- function (v) {
  p<- (-n*media/v)-(n)*log(v)-n*log(sum(exp(-x/v)))
  return(p) }

posterior4 <- function (v,phi) {
  p<--(n+1)*log(v)-(n*(media-phi)/v)-sum(exp(-(x-phi)/v))
  return(p)
}

posterior5 <- function (v,phi) {
  p <- dnorm(phi,mean=0,sd=1,log = TRUE)+
    dgamma(v,shape=4,rate=4,log = TRUE)-
    n*log(v)-(n*(media-phi)/v)-sum(exp(-(x-phi)/v))
  return(p)
}


fisher<-function(phi,var) {
  F<-matrix(nrow=2,ncol=2)
  F[1,1]<- 1/(var^2)
  F[1,2]<- -(1/(var^2))*(1+digamma(1))  #Meu
  F[2,1]<- F[1,2]
  F[2,2]<- (1/(var^2))*(1+trigamma(2)+digamma(2)^2)     #Meu
  return(F)
}

loglike<-function(theta) {
  phi<-theta[1]
  v<-theta[2]
  aux<--(n)*log(v)-(n*(media-phi)/v)-sum(exp(-(x-phi)/v))
  return(aux)
}


# Tamanho da amostra
n<-nini
for(k in 1:NN){
  emv <- matrix(nrow=B,ncol=2)
  emcmc1 <- matrix(nrow=B,ncol=2)
  emcmc2 <- matrix(nrow=B,ncol=2)
  emcmc3 <- matrix(nrow=B,ncol=2)
  emcmc4 <- matrix(nrow=B,ncol=2)
  pb11<-rep(0,times=B)
  pb12<-rep(0,times=B)
  pb21<-rep(0,times=B)
  pb22<-rep(0,times=B)
  pb31<-rep(0,times=B)
  pb32<-rep(0,times=B)
  pb41<-rep(0,times=B)
  pb42<-rep(0,times=B)
  ge1<-rep(0,times=B)
  ge2<-rep(0,times=B)
  ge3<-rep(0,times=B)
  ge4<-rep(0,times=B)
  ge5<-rep(0,times=B)
  ge6<-rep(0,times=B)
  eulerg = 0.57721566490153
  o<-0           # Contador de itera??es
  ite<-0         # Contador de itera??es
  
  ################## Calculando #######################
  ### Comando try() faz com que mesmo que o programa encontre erro ele n?o pare de executar ###
  while(o<B)
  {
    phi<-length(R+1)
    varp<-length(R+1)
    x<-rgumbel(n,pphi,pvarp)
    media<-mean(x)
    b<-1 #Valor Auxiliar MCMC
    m<-1
    iniphi<-media-(eulerg*sd(x)*sqrt(6)/pi)
    inivarp<-sd(x)*sqrt(6)/pi
    phi[1]<-iniphi  #Valores Iniciais
    varp[1]<-inivarp  #Valores Iniciais
    c1<-rep(0,times=R)    #Contador de valores aceitos
    c2<-rep(0,times=R)    ## Realizando o M-H hibrido
    a1<-0;  a2<-0;i<-1; c10<-0
    try(
      while (i<=R) {
        if(i<1) i<-2
        prop1<-rgamma(1,shape=b*varp[i],rate=b)
        ratio1<-posterior1(prop1)-posterior1(varp[i])+dgamma(varp[i],shape=b*prop1,rate=b,log=TRUE)-dgamma(prop1,shape=b*varp[i],rate=b,log=TRUE)
        has<-min(1,exp(ratio1)); u1<-runif(1)
        if (u1<has & is.double(has) & has!="NaN" & has!="Inf" & has!="-Inf" & prop1>0.01) {varp[i+1]<-prop1 ; c1[i]<-0 ; a1<-0} else {varp[i+1]<-varp[i] ; a1<-a1+1; c1[i]<-1}
        prop2<-rnorm(1,mean=phi[i],sd=m)
        ratio2<-posterior2(varp[i+1],prop2)-posterior2(varp[i+1],phi[i])+dnorm(phi[i],mean=prop2,sd=m,log=TRUE)-dnorm(prop2,mean=phi[i],sd=m,log=TRUE)
        alpha2<-min(1,exp(ratio2))
        u2<-runif(1)
        if (u2<alpha2 & alpha2!="NaN" & alpha2!="Inf" & alpha2!="-Inf") {phi[i+1]<-prop2 ;c2[i]<-0 ; a2<-0} else { phi[i+1]<-phi[i] ;
        a2<-a2+1; c2[i]<-1 }
        if(a1==50 | a2==50) {i<-i-50; a1=0 ; a2=0}
        i<-i+1 ; c10<-c10+1
        #if(c10==50000) {i=R+1; phi<-rep(0,times=R); varp<-rep(0,times=R)}
      }
    )
    try(vvarp<-varp[seq(burnin,R,jump)])
    try(vphi<- phi[seq(burnin,R,jump)])
    if(is.double(vvarp))
      if(sum(vvarp)!=0){
        old  =  options(digits=5)
        ace1<- (1-sum(c1)/length(c1))
        atc1<- mean(acf(vvarp,plot=F)$acf)
        atc2<- mean(acf(vphi,plot=F)$acf)
        ge1[o+1]<-abs(geweke.diag(vvarp)$z[1])
        ge2[o+1]<-abs(geweke.diag(vphi)$z[1])
        prai1<-quantile(vphi, probs = 0.025, na.rm = FALSE,names = FALSE,type = 7)
        pras1<-quantile(vphi, probs = 0.975, na.rm = FALSE,names = FALSE,type = 7)
        prli1<-quantile(vvarp, probs = 0.025, na.rm = FALSE,names = FALSE,type = 7)
        prls1<-quantile(vvarp, probs = 0.975, na.rm = FALSE,names = FALSE,type = 7)
        auxvarp1<-mean(vvarp); auxphi1<-mean(vphi); 
        if (ace1>=0.1 & ace1<=0.7 & auxphi1<20 & auxvarp1<20 & ge1[o+1]<1.96  & ge2[o+1]<1.96) {
          ############################Calculando a segunda posteriori ############################
          ########################################################################################
          phi<-length(R+1)
          varp<-length(R+1)
          media<-mean(x)
          phi[1]<-iniphi  #Valores Iniciais
          varp[1]<-inivarp  #Valores Iniciaiss
          c3<-rep(0,times=R)    #Contador de valores aceitos
          c4<-rep(0,times=R)
          ## Realizando o M-H hibrido
          a3<-0; a4<-0;i<-1; c10<-0
          try(
            while (i<=R) {
              if(i<1) i<-2
              prop3<-rgamma(1,shape=b*varp[i],rate=b)
              ratio3<-posterior3(prop3)-posterior3(varp[i])+dgamma(varp[i],shape=b*prop3,rate=b,log=TRUE)-dgamma(prop3,shape=b*varp[i],rate=b,log=TRUE)
              has<-min(1,exp(ratio3)); u3<-runif(1)
              if (u3<has & is.double(has) & has!="NaN" & has!="Inf" & has!="-Inf" & prop3>0.01) {varp[i+1]<-prop3 ; c3[i]<-0 ; a3<-0} else {varp[i+1]<-varp[i] ; a3<-a3+1; c3[i]<-1}
              
              prop4<-rnorm(1,mean=phi[i],sd=m)
              ratio4<-posterior4(varp[i+1],prop4)-posterior4(varp[i+1],phi[i])+dnorm(phi[i],mean=prop4,sd=m,log=TRUE)-dnorm(prop4,mean=phi[i],sd=m,log=TRUE)
              alpha4<-min(1,exp(ratio4))
              u4<-runif(1)
              
              if (u4<alpha4 & alpha4!="NaN" & alpha4!="Inf" & alpha4!="-Inf") {phi[i+1]<-prop4 ;c4[i]<-0 ; a4<-0} else { phi[i+1]<-phi[i] ;
              a4<-a4+1; c4[i]<-1 }
              if(a3==50 | a4==50) {i<-i-50; a3=0 ; a4=0}
              i<-i+1 ; c10<-c10+1
              if(c10==50000) {i=R+1; phi<-rep(0,times=R); varp<-rep(0,times=R)}
            }
          )
          vphi<-c()
          vvarp<-c()      
          try(vvarp<-varp[seq(burnin,R,jump)])
          try(vphi<- phi[seq(burnin,R,jump)])
          if(is.double(vvarp))
            if(sum(vvarp)!=0){
              old  =  options(digits=5)
              ace3<- (1-sum(c3)/length(c3))
              atc3<- mean(acf(vvarp,plot=F)$acf)
              atc4<- mean(acf(vphi,plot=F)$acf)
              ge3[o+1]<-abs(geweke.diag(vvarp)$z[1])
              ge4[o+1]<-abs(geweke.diag(vphi)$z[1])
              prai2<-quantile(vphi, probs = 0.025, na.rm = FALSE,names = FALSE,type = 7)
              pras2<-quantile(vphi, probs = 0.975, na.rm = FALSE,names = FALSE,type = 7)
              prli2<-quantile(vvarp, probs = 0.025, na.rm = FALSE,names = FALSE,type = 7)
              prls2<-quantile(vvarp, probs = 0.975, na.rm = FALSE,names = FALSE,type = 7)
              auxvarp2<-mean(vvarp); auxphi2<-mean(vphi); 
              if (ace3>=0.1 & ace3<=0.7 & auxphi2<20 & auxvarp2<20 & ge3[o+1]<1.96  & ge4[o+1]<1.96) {
              ############################Calculando a tercera posteriori ############################
              ########################################################################################
              phi<-length(R+1)
              varp<-length(R+1)
              media<-mean(x)
              phi[1]<-iniphi  #Valores Iniciais
              varp[1]<-inivarp  #Valores Iniciaiss
              c5<-rep(0,times=R)    #Contador de valores aceitos
              c6<-rep(0,times=R)
              ## Realizando o M-H hibrido
              a5<-0; a6<-0;i<-1; c10<-0
              try(
                while (i<=R) {
                  if(i<1) i<-2

                  prop5<-rgamma(1,shape=b*varp[i],rate=b)
                  ratio5<-posterior5(prop5,phi[i])-posterior5(varp[1],phi[i])+dgamma(varp[i],shape=b*prop5,rate=b,log=TRUE)-dgamma(prop5,shape=b*varp[i],rate=b,log=TRUE)
                  has<-min(1,exp(ratio5)); u5<-runif(1)
                  if (u5<has & is.double(has) & has!="NaN" & has!="Inf" & has!="-Inf" & prop5>0.01) {varp[i+1]<-prop5 ; c5[i]<-0 ; a5<-0} else {varp[i+1]<-varp[i] ; a5<-a5+1; c5[i]<-1}
                  
                  prop6<-rnorm(1,mean=phi[i],sd=m)
                  ratio6<-posterior5(varp[i+1],prop6)-posterior5(varp[i+1],phi[i])+dnorm(phi[i],mean=prop6,sd=m,log=TRUE)-dnorm(prop6,mean=phi[i],sd=m,log=TRUE)
                  alpha6<-min(1,exp(ratio6))
                  u6<-runif(1)
                  
                  if (u6<alpha6 & alpha6!="NaN" & alpha6!="Inf" & alpha6!="-Inf") {phi[i+1]<-prop6 ;c6[i]<-0 ; a6<-0} else { phi[i+1]<-phi[i]; a6<-a6+1; c6[i]<-1 } ;
                  
                  if(a5==50 | a6==50) {i<-i-50; a5=0 ; a6=0}
                  i<-i+1 ; c10<-c10+1
                  if(c10==50000) {i=R+1; phi<-rep(0,times=R); varp<-rep(0,times=R)}
                }
              )
              vphi<-c()
              vvarp<-c()      
              try(vvarp<-varp[seq(burnin,R,jump)])
              try(vphi<- phi[seq(burnin,R,jump)])
              if(is.double(vvarp))
                if(sum(vvarp)!=0){
                  old  =  options(digits=5)
                  ace3<- (1-sum(c3)/length(c3))
                  atc3<- mean(acf(vvarp,plot=F)$acf)
                  atc4<- mean(acf(vphi,plot=F)$acf)
                  ge5[o+1]<-abs(geweke.diag(vvarp)$z[1])
                  ge6[o+1]<-abs(geweke.diag(vphi)$z[1])
                  prai3<-quantile(vphi, probs = 0.025, na.rm = FALSE,names = FALSE,type = 7)
                  pras3<-quantile(vphi, probs = 0.975, na.rm = FALSE,names = FALSE,type = 7)
                  prli3<-quantile(vvarp, probs = 0.025, na.rm = FALSE,names = FALSE,type = 7)
                  prls3<-quantile(vvarp, probs = 0.975, na.rm = FALSE,names = FALSE,type = 7)
                  auxvarp3<-mean(vvarp); auxphi3<-mean(vphi); 
              if (ace3>=0.1 & ace3<=0.7 & auxphi3<20 & auxvarp3<20 & ge5[o+1]<1.96  & ge6[o+1]<1.96) {
                out<-NULL  
                res  <- try(maxNR(loglike, start=c(pphi, pvarp)))
                out[1]<-try(res$estimate[1])
                out[2]<-try(res$estimate[2])
                if(is.double(out[1]) & is.double(out[2])){
                  FM<- try(n*fisher(out[1],out[2]))
                  IVFM<-try(solve(FM))
                  vflam<-try(IVFM[1,1])
                  vfgama<-try(IVFM[2,2])
                  if(is.double(vflam))
                    if(is.double(out[1]) & vflam>0 & vfgama>0 & out[2]>0){
                      o<-(o+1); 
                      cc<-qnorm(0.975,0,1)
                      iic<-out[1]-cc*(vflam^0.5)  ; isc<-out[1]+cc*(vflam^0.5) ;   
                      iigama<-out[2]-cc*(vfgama^0.5)  ; isgama<-out[2]+cc*(vfgama^0.5) ;
                      emcmc1[o,1] <- auxphi1; emcmc1[o,2]<-auxvarp1;
                      emcmc2[o,1] <- auxphi2; emcmc2[o,2]<-auxvarp2;
                      emcmc3[o,1] <- auxphi3; emcmc3[o,2]<-auxvarp3;  
                      emcmc4[o,1] <- out[1]; emcmc4[o,2]<-out[2];  
                      {
                        if(prai1<=pphi & pras1>= pphi) pb11[o]<-1
                        if(prli1<=pvarp & prls1>= pvarp) pb12[o]<-1
                        if(prai2<=pphi & pras2>= pphi) pb21[o]<-1
                        if(prli2<=pvarp & prls2>= pvarp) pb22[o]<-1
                        if(prai3<=pphi & pras3>= pphi) pb31[o]<-1
                        if(prli3<=pvarp & prls3>= pvarp) pb32[o]<-1
                        if(iic<=pphi & isc>= pphi) pb41[o]<-1
                        if(iigama<=pvarp & isgama>= pvarp) pb42[o]<-1
                      };
                      cat(o,"     ",ite,"     ",n,"\n")
                      
                    }
                } }   } } } }
        }
    
    ite<-ite+1
    
  } 
  
  
  n<-n+jumpn
  
  
  
  ############# Imprimi #######################
  
  gematriz1 [k,]<-ge1
  gematriz2 [k,]<-ge2
  gematriz3 [k,]<-ge3
  gematriz4 [k,]<-ge4
  gematriz5 [k,]<-ge5
  gematriz6 [k,]<-ge6
  mediabayes1<-c((sum(emcmc1[,1]-pphi)/B),var(emcmc1[,1]),(sum(emcmc1[,2]-pvarp)/B),var(emcmc1[,2]))
  mediabayes2<-c((sum(emcmc2[,1]-pphi)/B),var(emcmc2[,1]),(sum(emcmc2[,2]-pvarp)/B),var(emcmc2[,2]))
  mediabayes3<-c((sum(emcmc3[,1]-pphi)/B),var(emcmc3[,1]),(sum(emcmc3[,2]-pvarp)/B),var(emcmc3[,2]))
  mediabayes4<-c((sum(emcmc4[,1]-pphi)/B),var(emcmc4[,1]),(sum(emcmc4[,2]-pvarp)/B),var(emcmc4[,2]))
  mrephi[k,]<-c(mediabayes1[1],mediabayes2[1],mediabayes3[1],mediabayes4[1])
  msephi[k,]<-c(mediabayes1[2],mediabayes2[2],mediabayes3[2],mediabayes4[2])
  cobphi[k,]<-c(round(sum(pb11)/B,3),round(sum(pb21)/B,3),round(sum(pb31)/B,3),round(sum(pb41)/B,3))
  
  mrevarp[k,]<-c(mediabayes1[3],mediabayes2[3],mediabayes3[3],mediabayes4[3])
  msevarp[k,]<-c(mediabayes1[4],mediabayes2[4],mediabayes3[4],mediabayes4[4])
  cobvarp[k,]<-c(round(sum(pb12)/B,3),round(sum(pb22)/B,3),round(sum(pb32)/B,3),round(sum(pb42)/B,3))
  
  ya<-seq(nini,nmax,jumpn)   
  
  par(mfrow=c(2,3))
  par(mai=c(0.6,0.6,0.1, 0.1)) 
  plot(ya,mrephi[,1],xlim=c(min(ya),max(ya)),ylim=c(-0.23,0.1),xlab="n",ylab=expression(paste("Bias"," (",phi,")")),
       type="l", col="dodgerblue3",lwd=2,lty=1,main=NULL)
  abline(h = 0, v = 0, col = "gray60")
  points(ya,mrephi[,2],type="l",col="goldenrod3",lwd=2,lty=1)
  points(ya,mrephi[,4],type="l",col="firebrick2",lwd=2,lty=1)
  points(ya,mrephi[,3],type="l",col="darkgreen",lwd=2,lty=1)
  
  
  par(mai=c(0.6,0.6,0.1, 0.1)) 
  plot(ya,msephi[,1],xlim=c(min(ya),max(ya)),ylim=c(0,1),xlab="n",ylab=expression(paste("MSE"," (",phi,")")),
       type="l", col="dodgerblue3",lwd=2,lty=1,main=NULL)
  abline(h = 0, v = 0, col = "gray60")
  points(ya,msephi[,4],type="l",col="firebrick2",lwd=2,lty=1)
  points(ya,msephi[,2],type="l",col="goldenrod3",lwd=2,lty=1)
  points(ya,msephi[,3],type="l",col="darkgreen",lwd=2,lty=1)
  
  legend(30,0.92,c("Jeffreys","Reference","Gamma","MLE"),lwd=c(2,2,2,2),lty=c(1,1,1,1),cex=1.4,col = c("dodgerblue4","goldenrod3","darkgreen","firebrick2"))
  
  par(mai=c(0.6,0.6,0.1, 0.1)) 
  plot(ya,cobphi[,1],xlim=c(min(ya),max(ya)),ylim=c(0.84,0.98),xlab="n",ylab=expression(paste("CP"," (",phi,")")),type="l", col="dodgerblue3",lwd=2,lty=1,main=NULL)
  abline(h = 0.95, v = 0.95, col = "gray60")
  points(ya,cobphi[,2],type="l",col="goldenrod3",lwd=2,lty=1)
  points(ya,cobphi[,4],type="l",col="firebrick2",lwd=2,lty=1)
  points(ya,cobphi[,3],type="l",col="darkgreen",lwd=2,lty=1)
  
  par(mai=c(0.6,0.6,0.1, 0.1)) 
  plot(ya,mrevarp[,1],xlim=c(min(ya),max(ya)),ylim=c(-0.5,0.5),xlab="n",ylab=expression(paste("Bias"," (",varphi,")")),type="l", col="dodgerblue3",lwd=2,lty=1,main=NULL)
  abline(h = 0, v = 0, col = "gray60")
  points(ya,mrevarp[,3],type="l",col="firebrick2",lwd=2,lty=1)
  points(ya,mrevarp[,2],type="l",col="goldenrod3",lwd=2,lty=1)
  points(ya,mrevarp[,4],type="l",col="darkgreen",lwd=2,lty=1)

  par(mai=c(0.6,0.6,0.1, 0.1)) 
  plot(ya,msevarp[,1],xlim=c(min(ya),max(ya)),ylim=c(0,1),xlab="n",ylab=expression(paste("MSE"," (",varphi,")")),type="l", col="dodgerblue3",lwd=2,lty=1,main=NULL)
  abline(h = 0, v = 0, col = "gray60")
  points(ya,msevarp[,4],type="l",col="firebrick2",lwd=2,lty=1)
  points(ya,msevarp[,2],type="l",col="goldenrod3",lwd=2,lty=1)
  points(ya,msevarp[,3],type="l",col="darkgreen",lwd=2,lty=1)
  
  par(mai=c(0.6,0.6,0.1, 0.1)) 
  plot(ya,cobvarp[,1],xlim=c(min(ya),max(ya)),ylim=c(0.76,1),xlab="n",ylab=expression(paste("CP"," (",varphi,")")),type="l", col="dodgerblue3",lwd=2,lty=1,main=NULL)
  abline(h = 0.95, v = 0.95, col = "gray60")
  points(ya,cobvarp[,2],type="l",col="goldenrod3",lwd=2,lty=1)
  points(ya,cobvarp[,4],type="l",col="firebrick2",lwd=2,lty=1)
  points(ya,cobvarp[,3],type="l",col="darkgreen",lwd=2,lty=1)
  
  
}

save.image("simu1.Rdata")

