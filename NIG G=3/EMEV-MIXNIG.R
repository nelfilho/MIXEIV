rm(list = ls())
options(mc.cores = parallel::detectCores())

#source("reiv.NIGmix.r")
source("/home/jeremias/Desktop/NIG - G=3/reiv.NIGmix.r") ## This line must be used for Ubuntu compilation

#source("FUNÇÕES NECESSÁRIAS.r")

# Libraries
library(rjags)
library(MASS)
library(mcmcplots)
library(mixsmsn)

#funçoes que precisam ser rodadas

#Arquivos dic ST
setwd("DIC_EMEmixST") #entrando na pasta DIC_EMEmixST
source("DICSTmixeiv.R")
source("dSTmixme.R")
source("dskewTmix.mult2.r")
source("dstni2.r")
setwd("..")

#Arquivos dic SSL
setwd("DIC_EMEmixSSL")
source("DICSSLmixeiv.r")
source("dskewSSL.r")
source("dskewSSLmix.mult.r")
source("dSSLmixme.r")

#Arquivos dic SN
setwd("..")
setwd("DIC_EMEmixSN")
source("DICSN.mixeiv.r")
source("dskewNmix.mult2.r")
source("dSNmixme.r")
source("dSNni2.r")

#Arquivos dic T
setwd("..")
setwd("DIC_EMEmixT")
source("DICT.mixeiv.r")
source("dTmix.mult2.r")
#source("dstni2.r")
source("dTmixme.r")

#Arquivos dic N
setwd("..")
setwd("DIC_EMEmixN")
source("DICNORmixeiv.r")
source("dnormix.mult.r")
source("dNORmixme.r")

#Arquivos dic NC
setwd("..")
setwd("DIC_EMEmixNC")
source("dCNA.r")
source("dskewNCmix.mult2.r")
source("dNCmixme.r")
source("DICNC.mixeiv.r")

#Arquivos dic SNC
setwd("..")
setwd("DIC_EMEmixSNC")
source("dsnci2.r")
source("dskewNCmix.mult2.r")
source("dSTmixme.r")
source("DICSNCmixeiv.r")
setwd("..")

#Arquivos dic SL
setwd("DIC_EMEmixSL")
source("DICSL.mixeiv.r")
source("dSLmixme.r")
source("dSLmix.mult2.r")
source("dSL.r")
setwd("..")
# Geração do modelo MEV-NIGmix-----------------------------------------------------------------
p  <- 3
G  <- 3
mu <- c(-10,1,10)
lambda<-c(-2,1,-2)
delta <- 0.5
gama <- 1 
pii <- c(0.3,0.3,0.4)
Omega <- diag(c(1,1,1)) 
alpha <- c(0.4,0.1)
beta <- c(0.8,0.9)
n <- 100

library(robustbase)    
y.teste <- reiv.NIGmix(n,alpha,beta,Omega,mu,lambda,gama,delta,pii)
Z=y.teste$Z
#----------------------------------------------------------------------------------------------
#valores iniciais
ini <- list()

peso=c(0.6,0.2,0.2)
a=c(NA,0.6, 0.3)
b=c(NA,0.4,0.7)
iOmega=c(8,3,1)
mu=c(2,8,-8)
Delta=c(0,4,2)
igamma2=c(8,6,1)
nu=4
w=sample(G,size=n,replace=TRUE,prob=peso)

#numeros de interções mcmc
niter <- 25000

#burning
nburn <- 5000

#lag
nthin <- 5

#numero de cadeis
n.chains <- 1

#Modelol NC ------------------------------------------------------------------------------------
dados <- list(z=Z,n=n, G=G, p=p,theta=rep(1,G))
ini <-list(peso=peso,a=a,b=b,iOmega=iOmega,mu=mu,igamma2=igamma2,w=w)
par.inf <- c("a","b","Omega","mu","gamma2","peso","rho","tau")

m.jags.NC <- jags.model("EMEmixNC.txt", data = dados, quiet = TRUE, n.chains = n.chains,inits=ini)
Scoda.NC  <- coda.samples(m.jags.NC, par.inf, n.iter = niter - nburn, thin=nthin)
datNC     <- as.data.frame(as.matrix(Scoda.NC))

nc_a=nc_b=nc_Omega=nc_mu=nc_gamma2=nc_peso=NULL
for(i in 1:p){
  nc_a[i] =  which(names(datNC)==paste("a[",i,"]",sep=''))
  nc_b[i] =  which(names(datNC)==paste("b[",i,"]",sep=''))
  nc_Omega[i] = which(names(datNC)==paste("Omega[",i,"]",sep=''))
}
for(i in 1:G){
  nc_mu[i] =  which(names(datNC)==paste("mu[",i,"]",sep=''))
  nc_gamma2[i] = which(names(datNC)==paste("gamma2[",i,"]",sep='')) 
  nc_peso[i] = which(names(datNC)==paste("peso[",i,"]",sep='')) 
}

MAT.NC = data.frame(parametros=rownames(summary(Scoda.NC)[[1]]),round(cbind(summary(Scoda.NC)$statistics,HPDinterval(Scoda.NC)[[1]][,1:2]),4))

dicNC = DICNC.mixeiv(Z,alpha=t(datNC[,nc_a[-1]]),beta=t(datNC[,nc_b[-1]]),mu=t(datNC[,nc_mu]),gama=t(datNC[,nc_gamma2]),
                     peso=t(datNC[,nc_peso]),Omega=t(datNC[,nc_Omega]),rho=datNC$"rho",tau=datNC$"tau")

#Modelol SNC -----------------------------------------------------------------------------------
dados <- list(z=Z,n=n, G=G, p=p,theta=rep(1,G))
ini <-list(peso=peso,a=a,b=b,iOmega=iOmega,mu=mu,Delta=Delta,igamma2=igamma2,w=w)
par.inf <- c("a","b","Omega","mu","Delta","gamma2","peso","rho","tau")

m.jags.SNC <- jags.model("EMEmixSNC.txt", data = dados, quiet = TRUE, n.chains = n.chains,inits=ini)
Scoda.SNC  <- coda.samples(m.jags.SNC, par.inf, n.iter = niter - nburn, thin=nthin)
datSNC     <- as.data.frame(as.matrix(Scoda.SNC))

snc_a=snc_b=snc_Omega=snc_mu=snc_gamma2=snc_peso=snc_Delta=NULL
for(i in 1:p){
  snc_a[i]     = which(names(datSNC)==paste("a[",i,"]",sep=''))
  snc_b[i]     = which(names(datSNC)==paste("b[",i,"]",sep=''))
  snc_Omega[i] = which(names(datSNC)==paste("Omega[",i,"]",sep=''))
}
for(i in 1:G){
  snc_mu[i]     = which(names(datSNC)==paste("mu[",i,"]",sep=''))
  snc_gamma2[i] = which(names(datSNC)==paste("gamma2[",i,"]",sep='')) 
  snc_peso[i]   = which(names(datSNC)==paste("peso[",i,"]",sep=''))
  snc_Delta[i]  = which(names(datSNC)==paste("Delta[",i,"]",sep=''))
}

MAT.SNC = data.frame(parametros=rownames(summary(Scoda.SNC)[[1]]),round(cbind(summary(Scoda.SNC)$statistics,HPDinterval(Scoda.SNC)[[1]][,1:2]),4))

dicSNC = DICSNC.mixeiv(Z,alpha=t(datSNC[,snc_a[-1]]),beta=t(datSNC[,snc_b[-1]]),mu=t(datSNC[,snc_mu]),gama=t(datSNC[,snc_gamma2]),
                       Delta=t(datSNC[,snc_Delta]),peso=t(datSNC[,snc_peso]),Omega=t(datSNC[,snc_Omega]),tau=datSNC$"tau",rho=datSNC$"rho")

#Modelol SSL ------------------------------------------------------------------------------------
dados <- list(z=Z,n=n, G=G, p=p,theta=rep(1,G))
ini <-list(peso=peso,a=a,b=b,iOmega=iOmega,mu=mu,Delta=Delta,igamma2=igamma2,nu=nu,w=w)
par.inf <- c("a","b","Omega","mu","Delta","gamma2","nu","peso")

m.jags.SSL <- jags.model("EMEmixSSL.txt", data = dados, quiet = TRUE, n.chains = n.chains,inits=ini)
Scoda.SSL <- coda.samples(m.jags.SSL, par.inf, n.iter = niter - nburn, thin=nthin)
datSSL <- as.data.frame(as.matrix(Scoda.SSL))

sll_a=sll_b=sll_Omega=sll_mu=sll_gamma2=sll_peso=sll_Delta=NULL
for(i in 1:p){
  sll_a[i]     = which(names(datSSL)==paste("a[",i,"]",sep=''))
  sll_b[i]     = which(names(datSSL)==paste("b[",i,"]",sep=''))
  sll_Omega[i] = which(names(datSSL)==paste("Omega[",i,"]",sep=''))
}
for(i in 1:G){
  sll_mu[i]     = which(names(datSSL)==paste("mu[",i,"]",sep=''))
  sll_gamma2[i] = which(names(datSSL)==paste("gamma2[",i,"]",sep='')) 
  sll_peso[i]   = which(names(datSSL)==paste("peso[",i,"]",sep=''))
  sll_Delta[i]  = which(names(datSSL)==paste("Delta[",i,"]",sep=''))
}

MAT.SSL = data.frame(parametros=rownames(summary(Scoda.SSL)[[1]]),round(cbind(summary(Scoda.SSL)$statistics,HPDinterval(Scoda.SSL)[[1]][,1:2]),4))

dicSSL = DICSSL.mixeiv(Z,alpha=t(datSSL[,sll_a[-1]]),beta=t(datSSL[,sll_b[-1]]),mu=t(datSSL[,sll_mu]),gama=t(datSSL[,sll_gamma2]),
                       Delta=t(datSSL[,sll_Delta]),peso=t(datSSL[,sll_peso]),Omega=t(datSSL[,sll_Omega]),nu=datSSL$"nu")

#Modelol SL ------------------------------------------------------------------------------------
dados <- list(z=Z,n=n, G=G, p=p,theta=rep(1,G))
ini <-list(peso=peso,a=a,b=b,iOmega=iOmega,mu=mu,igamma2=igamma2,nu=nu,w=w)
par.inf <- c("a","b","Omega","mu","gamma2","nu","peso")

m.jags.SL <- jags.model("EMEmixSL.txt", data = dados, quiet = TRUE, n.chains = n.chains,inits=ini)
Scoda.SL <- coda.samples(m.jags.SL, par.inf, n.iter = niter - nburn, thin=nthin)
datSL <- as.data.frame(as.matrix(Scoda.SL))

sl_a=sl_b=sl_Omega=sl_mu=sl_gamma2=sl_peso=NULL
for(i in 1:p){
  sl_a[i] =  which(names(datSL)==paste("a[",i,"]",sep=''))
  sl_b[i] =  which(names(datSL)==paste("b[",i,"]",sep=''))
  sl_Omega[i] = which(names(datSL)==paste("Omega[",i,"]",sep=''))
}
for(i in 1:G){
  sl_mu[i] =  which(names(datSL)==paste("mu[",i,"]",sep=''))
  sl_gamma2[i] = which(names(datSL)==paste("gamma2[",i,"]",sep='')) 
  sl_peso[i] = which(names(datSL)==paste("peso[",i,"]",sep='')) 
}

MAT.SL = data.frame(parametros=rownames(summary(Scoda.SL)[[1]]),round(cbind(summary(Scoda.SL)$statistics,HPDinterval(Scoda.SL)[[1]][,1:2]),4))

dicSL = DICSL.mixeiv(Z,alpha=t(datSL[,sl_a[-1]]),beta=t(datSL[,sl_b[-1]]),mu=t(datSL[,sl_mu]),gama=t(datSL[,sl_gamma2]),
                     peso=t(datSL[,sl_peso]),Omega=t(datSL[,sl_Omega]),nu=datSL$"nu")

#Modelol ST-------------------------------------------------------------------------------------
ini <-list(peso=peso,a=a,b=b,iOmega=iOmega,mu=mu,Delta=Delta,igamma2=igamma2,nu=nu,w=w)
par.inf <- c("a","b","Omega","mu","Delta","gamma2","nu","peso")
dados <- list(z=Z,n=n, G=G, ll=0.04, mm=0.5, p=p,theta=rep(1,G))

m.jags.ST <- jags.model("EMEmixST.txt", data = dados, quiet = TRUE, n.chains = n.chains,inits=ini)
Scoda.ST <- coda.samples(m.jags.ST, par.inf, n.iter = niter - nburn, thin=nthin)
datST <- as.data.frame(as.matrix(Scoda.ST))

st_a=st_b=st_Omega=st_mu=st_gamma2=st_peso=st_Delta=NULL
for(i in 1:p){
  st_a[i]     = which(names(datST)==paste("a[",i,"]",sep=''))
  st_b[i]     = which(names(datST)==paste("b[",i,"]",sep=''))
  st_Omega[i] = which(names(datST)==paste("Omega[",i,"]",sep=''))
}
for(i in 1:G){
  st_mu[i]     = which(names(datST)==paste("mu[",i,"]",sep=''))
  st_gamma2[i] = which(names(datST)==paste("gamma2[",i,"]",sep='')) 
  st_peso[i]   = which(names(datST)==paste("peso[",i,"]",sep=''))
  st_Delta[i]  = which(names(datST)==paste("Delta[",i,"]",sep=''))
}

MAT.ST = data.frame(parametros=rownames(summary(Scoda.ST)[[1]]),round(cbind(summary(Scoda.ST)$statistics,HPDinterval(Scoda.ST)[[1]][,1:2]),4))

dicST = DICST.mixeiv(Z,alpha=t(datST[,st_a[-1]]),beta=t(datST[,st_b[-1]]),mu=t(datST[,st_mu]),gama=t(datST[,st_gamma2]),
                     Delta=t(datST[,st_Delta]),peso=t(datST[,st_peso]),Omega=t(datST[,st_Omega]),nu=datST$"nu")

#Modelol SN ------------------------------------------------------------------------------------
dados <- list(z=Z,n=n, G=G, p=p,theta=rep(1,G))
ini <-list(peso=peso,a=a,b=b,iOmega=iOmega,mu=mu,Delta=Delta,igamma2=igamma2,w=w)
par.inf <- c("a","b","Omega","mu","Delta","gamma2","peso") #parametros que vc vai querer inferencia

m.jags.SN <- jags.model("EMEmixSN.txt", data = dados, quiet = TRUE, n.chains = n.chains,inits=ini)
Scoda.SN <- coda.samples(m.jags.SN, par.inf, n.iter = niter - nburn, thin=nthin)
datSN <- as.data.frame(as.matrix(Scoda.SN))

sn_a=sn_b=sn_Omega=sn_mu=sn_gamma2=sn_peso=sn_Delta=NULL
for(i in 1:p){
  sn_a[i]     = which(names(datSN)==paste("a[",i,"]",sep=''))
  sn_b[i]     = which(names(datSN)==paste("b[",i,"]",sep=''))
  sn_Omega[i] = which(names(datSN)==paste("Omega[",i,"]",sep=''))
}
for(i in 1:G){
  sn_mu[i]     = which(names(datSN)==paste("mu[",i,"]",sep=''))
  sn_gamma2[i] = which(names(datSN)==paste("gamma2[",i,"]",sep='')) 
  sn_peso[i]   = which(names(datSN)==paste("peso[",i,"]",sep=''))
  sn_Delta[i]  = which(names(datSN)==paste("Delta[",i,"]",sep=''))
}

MAT.SN = data.frame(parametros=rownames(summary(Scoda.SN)[[1]]),round(cbind(summary(Scoda.SN)$statistics,HPDinterval(Scoda.SN)[[1]][,1:2]),4))

dicSN = DICSN.mixeiv(Z,alpha=t(datSN[,sn_a[-1]]),beta=t(datSN[,sn_b[-1]]),mu=t(datSN[,sn_mu]),gama=t(datSN[,sn_gamma2]),
                     Delta=t(datSN[,sn_Delta]),peso=t(datSN[,sn_peso]),Omega=t(datSN[,sn_Omega]))

#Modelol T -------------------------------------------------------------------------------------
dados <- list(z=Z,n=n, G=G, ll=0.04, mm=0.5, p=p,theta=rep(1,G))
ini <-list(peso=peso,a=a,b=b,iOmega=iOmega,mu=mu,igamma2=igamma2,nu=nu,w=w)
par.inf <- c("a","b","Omega","mu","gamma2","nu","peso")

m.jags.T <- jags.model("EMEmixT.txt", data = dados, quiet = TRUE, n.chains = n.chains,inits=ini)
Scoda.T <- coda.samples(m.jags.T, par.inf, n.iter = niter - nburn, thin=nthin)
datT <- as.data.frame(as.matrix(Scoda.T))

t_a=t_b=t_Omega=t_mu=t_gamma2=t_peso=NULL
for(i in 1:p){
  t_a[i] =  which(names(datT)==paste("a[",i,"]",sep=''))
  t_b[i] =  which(names(datT)==paste("b[",i,"]",sep=''))
  t_Omega[i] = which(names(datT)==paste("Omega[",i,"]",sep=''))
}
for(i in 1:G){
  t_mu[i] =  which(names(datT)==paste("mu[",i,"]",sep=''))
  t_gamma2[i] = which(names(datT)==paste("gamma2[",i,"]",sep='')) 
  t_peso[i] = which(names(datT)==paste("peso[",i,"]",sep='')) 
}

MAT.T = data.frame(parametros=rownames(summary(Scoda.T)[[1]]),round(cbind(summary(Scoda.T)$statistics,HPDinterval(Scoda.T)[[1]][,1:2]),4))

dicT = DICT.mixeiv(Z,alpha=t(datT[,t_a[-1]]),beta=t(datT[,t_b[-1]]),mu=t(datT[,t_mu]),gama=t(datT[,t_gamma2]),
                   peso=t(datT[,t_peso]),Omega=t(datT[,t_Omega]),nu=datT$"nu")

#Modelol N -------------------------------------------------------------------------------------
dados <- list(z=Z,n=n, G=G, p=p,theta=rep(1,G))
ini <-list(peso=peso,a=a,b=b,iOmega=iOmega,mu=mu,igamma2=igamma2,w=w)
par.inf <- c("a","b","Omega","mu","gamma2","peso")

m.jags.N <- jags.model("EMEmixN.txt", data = dados, quiet = TRUE, n.chains = n.chains,inits=ini)
Scoda.N <- coda.samples(m.jags.N, par.inf, n.iter = niter - nburn, thin=nthin)
datN <- as.data.frame(as.matrix(Scoda.N))

n_a=n_b=n_Omega=n_mu=n_gamma2=n_peso=NULL
for(i in 1:p){
  n_a[i] =  which(names(datN)==paste("a[",i,"]",sep=''))
  n_b[i] =  which(names(datN)==paste("b[",i,"]",sep=''))
  n_Omega[i] = which(names(datN)==paste("Omega[",i,"]",sep=''))
}
for(i in 1:G){
  n_mu[i] =  which(names(datN)==paste("mu[",i,"]",sep=''))
  n_gamma2[i] = which(names(datN)==paste("gamma2[",i,"]",sep='')) 
  n_peso[i] = which(names(datN)==paste("peso[",i,"]",sep='')) 
}

MAT.N = data.frame(parametros=rownames(summary(Scoda.N)[[1]]),round(cbind(summary(Scoda.N)$statistics,HPDinterval(Scoda.N)[[1]][,1:2]),4))

dicN = DICNOR.mixeiv(Z,alpha=t(datN[,n_a[-1]]),beta=t(datN[,n_b[-1]]),mu=t(datN[,n_mu]),gama=t(datN[,n_gamma2]),
                     peso=t(datN[,n_peso]),Omega=t(datN[,n_Omega]))

#gravando dados---------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
dir.create(paste("Simulação",format(Sys.time(),"%d_%b_%Y às %Hh%Mm"),sep="_"))
setwd(paste("Simulação",format(Sys.time(),"%d_%b_%Y às %Hh%Mm"),sep="_"))


#tabelas
write.matrix(MAT.ST,"matriz de resultados ST.txt",sep="\t")
write.matrix(MAT.SSL,"matriz de resultados SSL.txt",sep="\t")
write.matrix(MAT.SN,"matriz de resultados SN.txt",sep="\t")
write.matrix(MAT.T,"matriz de resultados T.txt",sep="\t")
write.matrix(MAT.N,"matriz de resultados N.txt",sep="\t")
write.matrix(MAT.SNC,"matriz de resultados SNC.txt",sep="\t")
write.matrix(MAT.NC,"matriz de resultados NC.txt",sep="\t")
write.matrix(MAT.SL,"matriz de resultados SL.txt",sep="\t")

#graficos
aa = which(names(datN)=="a[1]" | names(datN)=="b[1]")
denplot(Scoda.N, parms = c(names(datN[-aa]), collapse = FALSE, auto.layout = TRUE))
savePlot("Density plot N",type="jpeg");graphics.off()
traplot(Scoda.N, parms = c(names(datN[-aa])))
savePlot("Trace plot N",type="jpeg");graphics.off()

aa = which(names(datT)=="a[1]" | names(datT)=="b[1]")
denplot(Scoda.T, parms = c(names(datT[-aa]), collapse = FALSE, auto.layout = TRUE))
savePlot("Density plot T",type="jpeg");graphics.off()
traplot(Scoda.T, parms = c(names(datT[-aa])))
savePlot("Trace plot T",type="jpeg");graphics.off()

aa = which(names(datSN)=="a[1]" | names(datSN)=="b[1]")
denplot(Scoda.SN, parms = c(names(datSN[-aa]), collapse = FALSE, auto.layout = TRUE))
savePlot("Density plot SN",type="jpeg");graphics.off()
traplot(Scoda.SN, parms = c(names(datSN[-aa])))
savePlot("Trace plot SN",type="jpeg");graphics.off()

aa = which(names(datSSL)=="a[1]" | names(datSSL)=="b[1]")
denplot(Scoda.SSL, parms = c(names(datSSL[-aa]), collapse = FALSE, auto.layout = TRUE))
savePlot("Density plot SSL",type="jpeg");graphics.off()
traplot(Scoda.SSL, parms = c(names(datSSL[-aa])))
savePlot("Trace plot SSL",type="jpeg");graphics.off()

aa = which(names(datST)=="a[1]" | names(datST)=="b[1]")
denplot(Scoda.ST, parms = c(names(datST[-aa]), collapse = FALSE, auto.layout = TRUE))
savePlot("Density plot ST",type="jpeg");graphics.off()
traplot(Scoda.ST, parms = c(names(datST[-aa])))
savePlot("Trace plot ST",type="jpeg");graphics.off()

aa = which(names(datNC)=="a[1]" | names(datNC)=="b[1]")
denplot(Scoda.NC, parms = c(names(datNC[-aa]), collapse = FALSE, auto.layout = TRUE))
savePlot("Density plot NC",type="jpeg");graphics.off()
traplot(Scoda.NC, parms = c(names(datNC[-aa])))
savePlot("Trace plot NC",type="jpeg");graphics.off()

aa = which(names(datSNC)=="a[1]" | names(datSNC)=="b[1]")
denplot(Scoda.SNC, parms = c(names(datSNC[-aa]), collapse = FALSE, auto.layout = TRUE))
savePlot("Density plot SNC",type="jpeg");graphics.off()
traplot(Scoda.SNC, parms = c(names(datSNC[-aa])))
savePlot("Trace plot SNC",type="jpeg");graphics.off()

aa = which(names(datSL)=="a[1]" | names(datSL)=="b[1]")
denplot(Scoda.SL, parms = c(names(datSL[-aa]), collapse = FALSE, auto.layout = TRUE))
savePlot("Density plot SL",type="jpeg");graphics.off()
traplot(Scoda.SL, parms = c(names(datSL[-aa])))
savePlot("Trace plot SL",type="jpeg");graphics.off()

#dics
Dic = matrix(as.matrix(cbind(dicN,dicT,dicNC,dicSN,dicST,dicSSL,dicSNC,dicSL)),3,8,dimnames=list(c("DIC.obs","tauD","Llog.hat"),c("N","T","NC","SN","ST","SSL","SNC","SL")))

m =as.data.frame(Dic[1,])
v=NULL
for(i in 1:8){
  v[i] = which(m==c(sort(m)[i]))
}
write.matrix(Dic[,v],"DIC.txt",sep="\t")

save(list = ls(),file = "rdata.RData")
setwd("..")

#library(gridExtra)
#p <- ggplot(as.data.frame(x=y.teste$x, w=factor(y.teste$w)), aes(factor(y.teste$w),y.teste$x))
#grafico_2 = p+geom_boxplot() + xlab("Reference Instruments")+ ylab("")
#grafico_1 = qplot(y.teste$x, geom="histogram",fill=I("white"),alpha=I(1),col=I("black"),ylab="Density",xlab="Reference Instrument")
#grid.arrange(grafico_1 , grafico_2,ncol=2)

#traplot(Scoda.ST, parms = c("a[2]","a[3]","b[2]","b[3]","Omega[1]","Omega[2]"),main = expression(alpha[1],alpha[2],beta[1],beta[2],omega[1],omega[2],omega[3]) ,col="black")
