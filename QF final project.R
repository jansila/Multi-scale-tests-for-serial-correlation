# Exam project Quantitative Finance II
setwd("~/Dropbox/My IES/QF II")
# Libraries
library(wavelets)
library(MASS)
library(mvtnorm)
library(TSA)
library(fGarch)
library(ggplot2)
library(reshape)
library(gplots)
#important constants
T<-2000
m<-4
max.m<-round(log2(T))
max.m

dev.off()
plot(modwt(wn, filter="haar", n.levels = 4), main="MODWT of White noise")
ploteps(wn)
wd<-arima.sim(n = 2000, model = list(ar = c(0.8897), order=c(1,0,0)))
plot(modwt(wd, filter="haar", n.levels = 4))
ploteps(wd)

ar1<-rar1(1000)
par(mfrow=c(3,7))
ploteps(ar1[,1])
ploteps(ar1[,2])
ploteps(ar1[,3])
ploteps(ar1[,4])
ploteps(ar1[,5])
ploteps(ar1[,6])
ploteps(ar1[,7])
ploteps(ar1[,8])
ploteps(ar1[,9])
ploteps(ar1[,10])
ploteps(ar1[,11])
ploteps(ar1[,12])
ploteps(ar1[,13])
ploteps(ar1[,14])
ploteps(ar1[,15])
ploteps(ar1[,16])
ploteps(ar1[,17])
ploteps(ar1[,18])
ploteps(ar1[,19])
ploteps(ar1[,20])
ploteps(ar1[,21])




# Funkce pro vypocet epsilon 1..4 a ploteps

ploteps<-function(series){
  i<-seq(1,4,1)
  trans<-modwt(series, filter="haar", n.levels=4)
  eps1<-sum(trans@W$W1^2)/sum(series^2)
  eps2<-sum(trans@W$W2^2)/sum(series^2)
  eps3<-sum(trans@W$W3^2)/sum(series^2)
  eps4<-sum(trans@W$W4^2)/sum(series^2)
  eps<-c(eps1,eps2,eps3,eps4)
  plot(eps,xlim=c(1,4),ylim=c(0,0.55),pch=20, col="red", lwd=2,main="Variance ratio theoretical and measured")
  lines(1/2^i, pch=4, lwd=2)
}

#Testovy statistiky
#according to Corollary 11


GS_4coef<-function(eps_coef){
  GS<-vector()
  names(GS)<-c("GS1", "GS2", "GS3","GS4")
  GS[1]<-c(sqrt(4*T)*(eps_coef[1]-1/2))  
  GS[2]<-c(sqrt(32*T/3)*(eps_coef[2]-1/4))
  GS[3]<-c(sqrt(256*T/15)*(eps_coef[3]-1/8))
  GS[4]<-c(sqrt(2048*T/71)*(eps_coef[4]-1/16))
  return(GS)
}




GS1<-function(series){
  T<-length(series)
  trans<-modwt(series, filter="haar", n.levels=1)
  gs1<-sqrt(4*T)*((sum(trans@W$W1^2)/sum(series^2))-1/2)
  rejection<-c(gs1>qnorm(0.975)|gs1<qnorm(0.025),gs1>qnorm(0.995)|gs1<qnorm(0.005))
  rejection}

rejGS1<-function(series){
  a<-colSums(t(apply(series, 2,GS1)))/5000
return(a*100)       
}


sGS1<-function(series){
  T<-length(series)
  trans<-modwt(series, filter="haar", n.levels=1)
  gs1<-sqrt(4*T)*((sum(trans@W$W1^2)/sum(series^2))-1/2)
  g<-c(gs1>qnorm(0.99)|gs1<qnorm(0.01))
g}


rejsGS1<-function(series){
  a<-colSums(t(apply(series, 2,sGS1)))/5000
  return(a*100)       
}


sGS2<-function(series){
  T<-length(series)
  trans<-modwt(series, filter="haar", n.levels=2)
  gs1<-(sqrt(32*T/3)*((sum(trans@W$W2^2)/sum(series^2))-1/4))
  g<-c(gs1>qnorm(0.99)|gs1<qnorm(0.01))
  g}


GS2<-function(series){
  T<-length(series)
  trans<-modwt(series, filter="haar", n.levels=2)
  gs1<-(sqrt(32*T/3)*((sum(trans@W$W2^2)/sum(series^2))-1/4))
  rejection<-c(gs1>qnorm(0.975)|gs1<qnorm(0.025),gs1>qnorm(0.995)|gs1<qnorm(0.005))
  rejection}          


rejGS2<-function(series){
  a<-colSums(t(apply(series, 2,GS2)))/5000
  return(a*100)       
}


GS1M<-function(series){
  T<-length(series)
  trans<-modwt(series, filter="haar", n.levels=1)
  eps1<-sum(trans@W$W1^2)/sum(series^2)
  eps_coef<-c(eps1)
  names(eps_coef)<-c("eps1")
  eps_print<-eps_coef[1:1]
  GS_temp<-vector()
  GS_temp[1]<-c(sqrt(4*T)*(eps_coef[1]-1/2))  
  names(GS_temp)<-c("GS1")
  GS<-GS_temp
  rejection<-c(GS[1]>qchisq(0.95,df=1), GS[1]>qchisq(0.99,df=1))  
  return(rejection)
}

GSM2<-function(series){
  T<-length(series)
  trans<-modwt(series, filter="haar", n.levels=2)
  eps1<-sum(trans@W$W1^2)/sum(series^2)
  eps2<-sum(trans@W$W2^2)/sum(series^2)
  eps_coef<-c(eps1,eps2)
  names(eps_coef)<-c("eps1","eps2")
  eps_print<-eps_coef[1:m]
  GS_temp<-vector()
  GS_temp[1]<-c(sqrt(4*T)*(eps_coef[1]-1/2))  
  GS_temp[2]<-c(sqrt(32*T/3)*(eps_coef[2]-1/4))
  A <- matrix(data = c(1,-1/sqrt(6),-1/sqrt(6),1),nrow = 2)
  GS<-as.matrix(t(GS_temp)) %*% A %*% t(t(as.matrix(GS_temp)))
  rejection<-c(GS[1]>qchisq(0.95,df=2), GS[1]>qchisq(0.99,df=2))  
  return(rejection)
}

GSM3<-function(series){
  T<-length(series)
  trans<-modwt(series, filter="haar", n.levels=3)
  eps1<-sum(trans@W$W1^2)/sum(series^2)
  eps2<-sum(trans@W$W2^2)/sum(series^2)
  eps3<-sum(trans@W$W3^2)/sum(series^2)
  eps_coef<-c(eps1,eps2,eps3)
  names(eps_coef)<-c("eps1","eps2","eps3")
  eps_print<-eps_coef[1:m]
  GS_temp<-vector()
  GS_temp[1]<-c(sqrt(4*T)*(eps_coef[1]-1/2))  
  GS_temp[2]<-c(sqrt(32*T/3)*(eps_coef[2]-1/4))
  GS_temp[3]<-c(sqrt(256*T/15)*(eps_coef[3]-1/8))
  names(GS_temp)<-c("GS1", "GS2", "GS3")
  A<-matrix(c(1,-1/sqrt(6),-5/sqrt(60),-1/sqrt(6),1,2/sqrt(360),-5/sqrt(60),2/sqrt(360),1), ncol=3,nrow=3)
  GS<-as.matrix(t(GS_temp)) %*% A %*% t(t(as.matrix(GS_temp)))
  rejection<-c(GS[1]>qchisq(0.95,df=3), GS[1]>qchisq(0.99,df=3))
  return(rejection)
}

GSM33<-function(series){
  T<-length(series)
  trans<-modwt(series, filter="haar", n.levels=3)
  eps1<-sum(trans@W$W1^2)/sum(series^2)
  eps2<-sum(trans@W$W2^2)/sum(series^2)
  eps3<-sum(trans@W$W3^2)/sum(series^2)
  eps_coef<-c(eps1,eps2,eps3)
  names(eps_coef)<-c("eps1","eps2","eps3")
  GS_temp<-vector()
  GS_temp[1]<-c(sqrt(4*T)*(eps_coef[1]-1/2))  
  GS_temp[2]<-c(sqrt(32*T/3)*(eps_coef[2]-1/4))
  GS_temp[3]<-c(sqrt(256*T/15)*(eps_coef[3]-1/8))
  names(GS_temp)<-c("GS1", "GS2", "GS3")
  A<-matrix(c(1,-1/sqrt(6),-5/sqrt(60),-1/sqrt(6),1,2/sqrt(360),-5/sqrt(60),2/sqrt(360),1), ncol=3,nrow=3)
  GS<-as.matrix(t(GS_temp)) %*% A %*% t(t(as.matrix(GS_temp)))
  rejection<-c(GS[1]>qchisq(0.99,df=3))
  return(rejection)
}




Qkstat<-function(series,k){
shit<-c(Box.test(x=series, type=c("Box-Pierce"), lag=k)$p.value<0.05, Box.test(x=series,type=c("Box-Pierce"), lag=k)$p.value<0.01)
shit}


#### Series simulations
# 5000 replications
# 100,300,1000,5000 obs
norm<-function(T){
  a<-matrix(0,ncol=5000,nrow = T) 
  for(i in 1:5000){
    a[,i]<-rnorm(T)
  }
  return(a)
}

norm100<-norm(100)
norm300<-norm(300)
norm1000<-norm(1000)
norm5000<-norm(5000)

garchn<-function(T){
a<-matrix(0,ncol=5000,nrow = T)
for(i in 1:5000){
  a[,i]<-garch.sim(alpha=c(0.001,0.05), beta=c(0.9),n=T)
}
return(a)
}
garchn100<-garchn(100)
garchn300<-garchn(300)
garchn1000<-garchn(1000)
garchn5000<-garchn(5000)


mix<-function(T){
  n<-T*5000
  s<-as.logical(rbinom(n, 1, 1/2))
  x<-vector(length=n)
  x[s]<-rnorm(n, mean=0,sd=1)[s]
  x[!s]<-rnorm(n, mean=0, sd=sqrt(1/2))[!s]
  matrix(x,nrow=T, ncol = 5000)
}
plot(density(mix100[,1]))

mix100<-mix(100)
mix300<-mix(300)
mix1000<-mix(1000)
mix5000<-mix(5000)

vart<-function(T){
s<-seq(from=1, to=50, by=50/T)
a<-rnorm(T,mean=0, sd=s)
a}


trend100<-replicate(5000, vart(100))
trend300<-replicate(5000, vart(300))
trend1000<-replicate(5000, vart(1000))
trend5000<-replicate(5000, vart(5000))
###########
# RAR(1), RAR(2)
###########
## AR(1) ##
q<-seq(from=-0.5, to=0.5, by=0.05)


rar1<-function(T=1000){
  rar<-matrix(0, ncol = length(q), nrow = T)
 for(i in 1:length(q)){ 
rar[,i]<-arima.sim(n=T, list(ar=q[i]))
}
rar}

rar2<-function(T=1000){
  rar<-matrix(0, ncol = length(q), nrow = T)
  for(i in 1:length(q)){ 
    rar[,i]<-arima.sim(n=T, list(ar=c(0,q[i])))
  }
  rar}

rejsGS1<-function(T){
a<-matrix(0, ncol=length(q), nrow=1000)
for(i in 1:1000)
{
    rar<-rar1(T)
a[i,]<-apply(rar,2,sGS1)
}
colMeans(a)
}

rejsGS1.ar2<-function(T){
  a<-matrix(0, ncol=length(q), nrow=1000)
  for(i in 1:1000)
  {
    rar<-rar2(T)
    a[i,]<-apply(rar,2,sGS1)
  }
  colMeans(a)
}

rejsGS2<-function(T){
  a<-matrix(0, ncol=length(q), nrow=1000)
  for(i in 1:1000)
  {
    rar<-rar1(T)
    a[i,]<-apply(rar,2,sGS2)
  }
  colMeans(a)
}


rejsGS2.ar2<-function(T){
  a<-matrix(0, ncol=length(q), nrow=1000)
  for(i in 1:1000)
  {
    rar<-rar2(T)
    a[i,]<-apply(rar,2,sGS2)
  }
  colMeans(a)
}
rejGGSM3<-function(T){
  a<-matrix(0, ncol=length(q), nrow=1000)
  for(i in 1:1000)
  {
    rar<-rar1(T)
    a[i,]<-apply(rar,2,GSM33)
  }
  colMeans(a)
}

rejGGSM3.ar2<-function(T){
  a<-matrix(0, ncol=length(q), nrow=1000)
  for(i in 1:1000)
  {
    rar<-rar2(T)
    a[i,]<-apply(rar,2,GSM33)
  }
  colMeans(a)
}


t1<-rejsGS1(100)
t2<-rejsGS1(300)
t3<-rejsGS1(1000)

tt1<-rejsGS2(100)
tt2<-rejsGS2(300)
tt3<-rejsGS2(1000)

ttt1<-rejGGSM3(100)
ttt2<-rejGGSM3(300)
ttt3<-rejGGSM3(1000)



at1<-rejsGS1.ar2(100)
at2<-rejsGS1.ar2(300)
at3<-rejsGS1.ar2(1000)

att1<-rejsGS2.ar2(100)
att2<-rejsGS2.ar2(300)
att3<-rejsGS2.ar2(1000)

attt1<-rejGGSM3.ar2(100)
attt2<-rejGGSM3.ar2(300)
attt3<-rejGGSM3.ar2(1000)


#### FINAL DATA #####
ar1.gs1<-rbind(t1,t2,t3) #AR1 process GS statistics
ar1.gs2<-rbind(tt1,tt2,tt3)
ar1.gsm3<-rbind(tt1,tt2,tt3)

ar2.gs1<-rbind(at1,at2,at3) #AR2 process GS statistics
ar2.gs2<-rbind(att1,att2,att3)
ar2.gsm3<-rbind(att1,att2,att3)

#### FINAL GRAPH ########
dev.off()
par(mfrow=c(3,2))

###
###
###  Pridelat legendu!
###
###
plot(t1, x=q, ylim=c(0,1), pch=0, ylab="Power",xlab="AR(1) coefficient", main="GS1 test power", col="darkblue")
points(t1, pch=22, col="darkblue")
points(t2, x=q, pch=2,, col="darkred")
lines(t2, col="darkred")
points(t3, x=q, pch=22)
lines(t3)

plot(at1, x=q, ylim=c(0,1),pch=0,ylab="Power",xlab="AR(2) coefficient", main="GS1 test power", col="darkblue")
lines(at1, col="darkblue")
points(at2, x=q, pch=2, col="darkred")
lines(at2, col="darkred")
points(at3, x=q, pch=22)
lines(at3)

plot(tt1, x=q,ylim=c(0,1), pch=0,ylab="Power",xlab="AR(1) coefficient", main="GS2 test power", col="darkblue")
lines(tt1, col="darkblue")
points(tt2, x=q, pch=2, col="darkred")
lines(tt1, col="darkred")
points(tt3, x=q, pch=22)
lines(tt3)


plot(att1, x=q, ylim=c(0,1),pch=0,ylab="Power",xlab="AR(2) coefficient", main="GS2 test power", col="darkblue")
lines(att1, col="darkblue")
points(att2, x=q, pch=2, col="darkred")
lines(att2, col="darkred")
points(att3, x=q, pch=22)
lines(att3)

plot(ttt1, x=q, ylim=c(0,1),pch=0,ylab="Power", xlab="AR(1) coefficient", main="GSM3 test power", col="darkblue")
lines(ttt1, col="darkblue")
points(ttt2, x=q, pch=2, col="darkred")
lines(ttt2, col="darkred")
points(ttt3, x=q, pch=22)
lines(ttt3)


plot(attt1, x=q,ylim=c(0,1), pch=0,ylab="Power", xlab="AR(2) coefficient", main="GSM3 test power", col="darkblue")
lines(attt1, col="darkblue")
points(attt2, x=q, pch=2, col="darkred")
lines(attt2, col="darkred")
points(attt3, x=q, pch=22)
lines(attt3)

########## Final GRAPH end #######
####### FINAL TABLE  ##########

vysledky<-rbind(vysledkyGS1,vysledkyGS2,vysledkyGSM2,vysledkyGSM3,Qkk)
vysledky
write.csv2(vysledky, file="vysledky.csv")
save(vysledky,file="vysledky.csv")
?save
##### FINAL TABLE GRAPH end #######
#
# RUN THE CODE DOWN IF NOT BEFORE
#
#




###  pokus heat map
###
v<-melt(vysledky)
nba_heatmap <- heatmap(t(as.matrix(vysledky)), Rowv=NA, Colv=NA, col = heat.colors(122), scale="column", margins=c(5,10))
######### pokus heatmap ale prcat to


# REJRATES for GSM
rejGSM2<-function(series){
  a<-colSums(t(apply(series, 2,GSM2)))/5000
return(a*100)
}
rejGSM1<-function(series){
  a<-colSums(t(apply(series, 2,GSM1)))/5000
  return(a*100)
}
rejGSM3<-function(series){
  a<-colSums(t(apply(series, 2,GSM3)))/5000
  return(a*100)
}


Qkrej<-function(series,k){
rej<-matrix(0, ncol = 2, nrow = 5000)
for (i in 1:5000){
rej[i,]<-cbind(Box.test(x=series[,i], type=c("Box-Pierce"), lag=k)$p.value<0.05,Box.test(x=series[,i], type=c("Box-Pierce"), lag=k)$p.value<0.01)
}
apply(rej,2,mean)*100}
  
Qk<-function(k){
    tab51<-cbind(Qkrej(norm100,k),Qkrej(norm300,k),Qkrej(norm1000,k),Qkrej(norm5000,k))
    tab52<-cbind(Qkrej(garchn100,k),Qkrej(garchn300,k),Qkrej(garchn1000,k),Qkrej(garchn5000,k))
    tab53<-cbind(Qkrej(mix100,k),Qkrej(mix300,k),Qkrej(mix1000,k),Qkrej(mix5000,k))
    tab54<-cbind(Qkrej(trend100,k),Qkrej(trend300,k),Qkrej(trend1000,k),Qkrej(trend5000,k))
    
    vysledkyQ5<-matrix(t(c(tab51,tab52,tab53,tab54)), nrow = 2)
    colnames(vysledkyQ5)<-c("norm100","norm300","norm1000","norm5000","garchn100","garchn300","garchn1000","garchn5000","mix100","mix300","mix1000","mix5000","trend100","trend300","trend1000","trend5000")
    vysledkyQ5}



###########
# rejection rates Box Pierce


Qkk<-rbind(Qk(5),Qk(10),Qk(20))
colnames(Qkk)<-c("norm100","norm300","norm1000","norm5000","garchn100","garchn300","garchn1000","garchn5000","mix100","mix300","mix1000","mix5000","trend100","trend300","trend1000","trend5000")

rownames(Qkk)<-c("Q5:95conf.l","Q5:99conf.l","Q10:95conf.l","Q10:99conf.l","Q20:95conf.l","Q20:99conf.l")

### rejection GS1 
tab1<-cbind(rejGS1(norm100),rejGS1(norm300),rejGS1(norm1000),rejGS1(norm5000))
tab2<-cbind(rejGS1(garchn100),rejGS1(garchn300),rejGS1(garchn1000),rejGS1(garchn5000))
tab3<-cbind(rejGS1(mix100),rejGS1(mix300),rejGS1(mix1000),rejGS1(mix5000))
tab4<-cbind(rejGS1(trend100),rejGS1(trend300),rejGS1(trend1000),rejGS1(trend5000))

vysledkyGS1<-matrix(t(c(tab1,tab2,tab3,tab4)), nrow = 2)
colnames(vysledkyGS1)<-c("norm100","norm300","norm1000","norm5000","garchn100","garchn300","garchn1000","garchn5000","mix100","mix300","mix1000","mix5000","trend100","trend300","trend1000","trend5000")
rownames(vysledkyGS1)<-c("GS1:95conf.l","GS1:99conf.l")
## konec GS1

## rejection GS 2 
tab1<-cbind(rejGS2(norm100),rejGS2(norm300),rejGS2(norm1000),rejGS2(norm5000))
tab2<-cbind(rejGS2(garchn100),rejGS2(garchn300),rejGS2(garchn1000),rejGS2(garchn5000))
tab3<-cbind(rejGS2(mix100),rejGS2(mix300),rejGS2(mix1000),rejGS2(mix5000))
tab4<-cbind(rejGS2(trend100),rejGS2(trend300),rejGS2(trend1000),rejGS2(trend5000))

vysledkyGS2<-matrix(t(c(tab1,tab2,tab3,tab4)), nrow = 2)
colnames(vysledkyGS2)<-c("norm100","norm300","norm1000","norm5000","garchn100","garchn300","garchn1000","garchn5000","mix100","mix300","mix1000","mix5000","trend100","trend300","trend1000","trend5000")
rownames(vysledkyGS2)<-c("GS2:95conf.l","GS2:99conf.l")

## konec GS2
##  rejection rates GSM2

tab1<-cbind(rejGSM2(norm100),rejGSM2(norm300),rejGSM2(norm1000),rejGSM2(norm5000))
tab2<-cbind(rejGSM2(garchn100),rejGSM2(garchn300),rejGSM2(garchn1000),rejGSM2(garchn5000))
tab3<-cbind(rejGSM2(mix100),rejGSM2(mix300),rejGSM2(mix1000),rejGSM2(mix5000))
tab4<-cbind(rejGSM2(trend100),rejGSM2(trend300),rejGSM2(trend1000),rejGSM2(trend5000))

vysledkyGSM2<-matrix(t(c(tab1,tab2,tab3,tab4)), nrow = 2)
colnames(vysledkyGSM2)<-c("norm100","norm300","norm1000","norm5000","garchn100","garchn300","garchn1000","garchn5000","mix100","mix300","mix1000","mix5000","trend100","trend300","trend1000","trend5000")
rownames(vysledkyGSM2)<-c("GSM2:95conf.l","GSM2:99conf.l")

####### KONEC GSM2 ######

#### GSM3 rejections

tab1<-cbind(rejGSM3(norm100),rejGSM3(norm300),rejGSM3(norm1000),rejGSM3(norm5000))
tab2<-cbind(rejGSM3(garchn100),rejGSM3(garchn300),rejGSM3(garchn1000),rejGSM3(garchn5000))
tab3<-cbind(rejGSM3(mix100),rejGSM3(mix300),rejGSM3(mix1000),rejGSM3(mix5000))
tab4<-cbind(rejGSM3(trend100),rejGSM3(trend300),rejGSM3(trend1000),rejGSM3(trend5000))

vysledkyGSM3<-matrix(t(c(tab1,tab2,tab3,tab4)), nrow = 2)
colnames(vysledkyGSM3)<-c("norm100","norm300","norm1000","norm5000","garchn100","garchn300","garchn1000","garchn5000","mix100","mix300","mix1000","mix5000","trend100","trend300","trend1000","trend5000")
rownames(vysledkyGSM3)<-c("GSM3:95conf.l","GSM3:99conf.l")


###### KONEC GSM3 #####




###########################################################################################################
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
##########################################     END OF CORE CODE     #######################################
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
###########################################################################################################






i<-seq(1,10,1)
plot(1/2^(i))
points(1/2^(i+1))
