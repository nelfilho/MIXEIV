rm(list = ls())
getwd()

#library(rstan)
#library(mnormt)
library(msm)
library(MCMCpack)

data= read.table("lupos.txt",h=T)
Z = matrix(c(data$ptn24h_exame1,data$rpc_exame1),dim(data),dimnames=list(1:nrow(data),c("X","Y")))
Z = cbind((Z[,1])/1000,(Z[,2])/1000)
n = nrow(Z)

n.iter <- 25000
n.burn <- 5000
n.thin <- 5

source("funções do pacote mnormt.r")

dir.create(paste("Simulação",format(Sys.time(),"%d_%b_%Y às %Hh%Mm"),sep="_"))
setwd(paste("Simulação",format(Sys.time(),"%d_%b_%Y às %Hh%Mm"),sep="_"))

#Normal case

source("C:/Users/Nelson Lima/Dropbox/Artigo MIXEIV/Simulação e Banco de dados Stan e Jags/Lupos stan/DIC_EMEmixN/DICNORmixeiv.r")
source("C:/Users/Nelson Lima/Dropbox/Artigo MIXEIV/Simulação e Banco de dados Stan e Jags/Lupos stan/DIC_EMEmixN/dnormix.mult.r")
source("C:/Users/Nelson Lima/Dropbox/Artigo MIXEIV/Simulação e Banco de dados Stan e Jags/Lupos stan/DIC_EMEmixN/dNORmixme.r")
source("C:/Users/Nelson Lima/Dropbox/Artigo MIXEIV/Simulação e Banco de dados Stan e Jags/Lupos stan/estimate_N/Stan Model simple N.r")


fitN   <- rstan::stan(model_code=stan_codeN, iter=n.iter,warmup = n.burn, thin = n.thin, chains=1, 
              data = list(N = nrow(Z), xobs = Z[,1], y = Z[,2]),save_dso = FALSE)
nparN = 9
codaN = as.data.frame(fitN)
namesN = names(codaN)[1:nparN]
windows(width=1e+8, height=75e+5)
par(mfrow=c(3,3))
plot.ts(codaN$beta,main="beta")
plot.ts(codaN$alpha,main="alpha")
plot.ts(codaN$omega2_1,main="omega2_1")
plot.ts(codaN$omega2_0,main="omega2_0")
plot.ts(codaN$"mu[1]",main="mu[1]")
plot.ts(codaN$"mu[2]",main="mu[2]")
plot.ts(codaN$"gamma2[1]",main="gamma2[1]")
plot.ts(codaN$"gamma2[2]",main="gamma2[2]")
plot.ts(codaN$peso,main="peso")


savePlot("Traceplot modelo EMEmixN",type="jpeg");graphics.off()


dicN = DICNOR.mixeiv(Z,alpha=t(codaN$"alpha"),beta=t(codaN$"beta"),mu=t(cbind(codaN$"mu[1]",codaN$"mu[2]")),
                     gama=t(cbind(codaN$"gamma2[1]",codaN$"gamma2[2]")),peso=t(cbind(codaN$"peso",1-codaN$"peso")),
                     Omega=t(cbind(codaN$"omega2_0",codaN$"omega2_1")))

matN = matrix(c(namesN,round(cbind(colMeans(codaN),apply(codaN,2,"median"),apply(codaN,2,"sd"),HPDinterval(mcmc(codaN))),4)[1:nparN,]),nparN,6,
       dimnames=list(1:nparN,c("parameters","mean","median","sd","lower","upper")))

capture.output(as.data.frame(matN),file = "matriz de estimação- N.txt")
save(list = ls(),file = "rdata.RData")

# t case

source("C:/Users/Nelson Lima/Dropbox/Artigo MIXEIV/Simulação e Banco de dados Stan e Jags/Lupos stan/DIC_EMEmixT/DICT.mixeiv.r")
source("C:/Users/Nelson Lima/Dropbox/Artigo MIXEIV/Simulação e Banco de dados Stan e Jags/Lupos stan/DIC_EMEmixT/dTmix.mult2.r")
source("C:/Users/Nelson Lima/Dropbox/Artigo MIXEIV/Simulação e Banco de dados Stan e Jags/Lupos stan/DIC_EMEmixT/dTmixme.r")
source("C:/Users/Nelson Lima/Dropbox/Artigo MIXEIV/Simulação e Banco de dados Stan e Jags/Lupos stan/DIC_EMEmixT/dstni2.r")
source("C:/Users/Nelson Lima/Dropbox/Artigo MIXEIV/Simulação e Banco de dados Stan e Jags/Lupos stan/estimate_T/Stan Model simple T.r")


fitT   <- rstan::stan(model_code=stan_codeT, iter=n.iter,warmup = n.burn, thin = n.thin, chains=1, 
              data = list(N = nrow(Z), xobs = Z[,1], y = Z[,2]),save_dso = FALSE)
nparT = 10
codaT = as.data.frame(fitT)
namesT = names(codaT)[1:nparT]
windows(width=1e+8, height=75e+5)
par(mfrow=c(2,5))
plot.ts(codaT$beta,main="beta")
plot.ts(codaT$alpha,main="alpha")
plot.ts(codaT$omega2_1,main="omega2_1")
plot.ts(codaT$omega2_0,main="omega2_0")
plot.ts(codaT$"mu[1]",main="mu[1]")
plot.ts(codaT$"mu[2]",main="mu[2]")
plot.ts(codaT$"gamma2[1]",main="gamma2[1]")
plot.ts(codaT$"gamma2[2]",main="gamma2[2]")
plot.ts(codaT$peso,main="peso")
plot.ts(codaT$nu,main="nu") 

savePlot("Traceplot modelo EMEmixT",type="jpeg");graphics.off()

dicT = DICT.mixeiv(Z,alpha=t(codaT$"alpha"),beta=t(codaT$"beta"),mu=t(cbind(codaT$"mu[1]",codaT$"mu[2]")),
                   gama=t(cbind(codaT$"gamma2[1]",codaT$"gamma2[2]")),peso=t(cbind(codaT$"peso",1-codaT$"peso")),
                   Omega=t(cbind(codaT$"omega2_0",codaT$"omega2_1")),nu=codaT$"nu")

matT = matrix(c(namesT,round(cbind(colMeans(codaT),apply(codaT,2,"median"),apply(codaT,2,"sd"),HPDinterval(mcmc(codaT))),4)[1:nparT,]),nparT,6,
              dimnames=list(1:nparT,c("parameters","mean","median","sd","lower","upper")))

capture.output(as.data.frame(matT),file = "matriz de estimação- T.txt")
save(list = ls(),file = "rdata.RData")

# SL case

source("C:/Users/Nelson Lima/Dropbox/Artigo MIXEIV/Simulação e Banco de dados Stan e Jags/Lupos stan/DIC_EMEmixSL/DICSL.mixeiv.r")
source("C:/Users/Nelson Lima/Dropbox/Artigo MIXEIV/Simulação e Banco de dados Stan e Jags/Lupos stan/DIC_EMEmixSL/dSL.r")
source("C:/Users/Nelson Lima/Dropbox/Artigo MIXEIV/Simulação e Banco de dados Stan e Jags/Lupos stan/DIC_EMEmixSL/dSLmix.mult2.r")
source("C:/Users/Nelson Lima/Dropbox/Artigo MIXEIV/Simulação e Banco de dados Stan e Jags/Lupos stan/DIC_EMEmixSL/dSLmixme.r")
source("C:/Users/Nelson Lima/Dropbox/Artigo MIXEIV/Simulação e Banco de dados Stan e Jags/Lupos stan/estimate_SL/Stan Model simple SL.r")

initf = function(){list(nu=2, gamma2=c(7,0.1))}
fitSL   <- rstan::stan(model_code=stan_codeSL, iter=n.iter,warmup = n.burn, thin = n.thin, chains=1, 
               data = list(N = nrow(Z), xobs = Z[,1], y = Z[,2]),save_dso = FALSE,init=initf)
nparSL = 10
codaSL = as.data.frame(fitSL)
namesSL = names(codaSL)[1:nparSL]
windows(width=1e+8, height=75e+5)

par(mfrow=c(2,5))
plot.ts(codaSL$beta,main="beta")
plot.ts(codaSL$alpha,main="alpha")
plot.ts(codaSL$omega2_1,main="omega2_1")
plot.ts(codaSL$omega2_0,main="omega2_0")
plot.ts(codaSL$"mu[1]",main="mu[1]")
plot.ts(codaSL$"mu[2]",main="mu[2]")
plot.ts(codaSL$"gamma2[1]",main="gamma2[1]")
plot.ts(codaSL$"gamma2[2]",main="gamma2[2]")
plot.ts(codaSL$peso,main="peso")
plot.ts(codaSL$nu,main="nu") 

savePlot("Traceplot modelo EMEmixSL",type="jpeg");graphics.off()

dicSL = DICSL.mixeiv(Z,alpha=t(codaSL$"alpha"),beta=t(codaSL$"beta"),mu=t(cbind(codaSL$"mu[1]",codaSL$"mu[2]")),
                     gama=t(cbind(codaSL$"gamma2[1]",codaSL$"gamma2[2]")),peso=t(cbind(codaSL$"peso",1-codaSL$"peso")),
                     Omega=t(cbind(codaSL$"omega2_0",codaSL$"omega2_1")),nu=codaSL$"nu")

matSL = matrix(c(namesSL,round(cbind(colMeans(codaSL),apply(codaSL,2,"median"),apply(codaSL,2,"sd"),HPDinterval(mcmc(codaSL))),4)[1:nparSL,]),nparSL,6,
              dimnames=list(1:nparSL,c("parameters","mean","median","sd","lower","upper")))

capture.output(as.data.frame(matSL),file = "matriz de estimação- SL.txt")
save(list = ls(),file = "rdata.RData")

# CN case

source("C:/Users/Nelson Lima/Dropbox/Artigo MIXEIV/Simulação e Banco de dados Stan e Jags/Lupos stan/DIC_EMEmixCN/dCNA.r")
source("C:/Users/Nelson Lima/Dropbox/Artigo MIXEIV/Simulação e Banco de dados Stan e Jags/Lupos stan/DIC_EMEmixCN/DICNC.mixeiv.r")
source("C:/Users/Nelson Lima/Dropbox/Artigo MIXEIV/Simulação e Banco de dados Stan e Jags/Lupos stan/DIC_EMEmixCN/dNCmixme.r")
source("C:/Users/Nelson Lima/Dropbox/Artigo MIXEIV/Simulação e Banco de dados Stan e Jags/Lupos stan/DIC_EMEmixCN/dskewNCmix.mult2.r")
source("C:/Users/Nelson Lima/Dropbox/Artigo MIXEIV/Simulação e Banco de dados Stan e Jags/Lupos stan/estimate_CN/Stan Model simple CN2.r")


fitCN   <- rstan::stan(model_code=stan_codeCN2, iter=n.iter,warmup = n.burn, thin = n.thin, chains=1, 
                data = list(N = nrow(Z), xobs = Z[,1], y = Z[,2]),save_dso = FALSE)
nparCN = 11
codaCN = as.data.frame(fitCN)
namesCN = names(codaCN)[1:nparCN]
windows(width=1e+8, height=75e+5)

par(mfrow=c(2,6))
plot.ts(codaCN$beta,main="beta")
plot.ts(codaCN$alpha,main="alpha")
plot.ts(codaCN$omega2_1,main="omega2_1")
plot.ts(codaCN$omega2_0,main="omega2_0")
plot.ts(codaCN$"mu[1]",main="mu[1]")
plot.ts(codaCN$"mu[2]",main="mu[2]")
plot.ts(codaCN$"gamma2[1]",main="gamma2[1]")
plot.ts(codaCN$"gamma2[2]",main="gamma2[2]")
plot.ts(codaCN$peso,main="peso")
plot.ts(codaCN$tau,main="tau") 
plot.ts(codaCN$rho,main="rho") 

savePlot("Traceplot modelo EMEmixCN",type="jpeg");graphics.off()

dicCN = DICNC.mixeiv(Z,alpha=t(codaCN$"alpha"),beta=t(codaCN$"beta"),mu=t(cbind(codaCN$"mu[1]",codaCN$"mu[2]")),
                     gama=t(cbind(codaCN$"gamma2[1]",codaCN$"gamma2[2]")),peso=t(cbind(codaCN$"peso",1-codaCN$"peso")),
                     Omega=t(cbind(codaCN$"omega2_0",codaCN$"omega2_1")),rho=codaCN$"rho",tau=codaCN$"tau")

matCN = matrix(c(namesCN,round(cbind(colMeans(codaCN),apply(codaCN,2,"median"),apply(codaCN,2,"sd"),HPDinterval(mcmc(codaCN))),4)[1:nparCN,]),nparCN,6,
               dimnames=list(1:nparCN,c("parameters","mean","median","sd","lower","upper")))

capture.output(as.data.frame(matCN),file = "matriz de estimação- CN.txt")
save(list = ls(),file = "rdata.RData")

# SN case

source("C:/Users/Nelson Lima/Dropbox/Artigo MIXEIV/Simulação e Banco de dados Stan e Jags/Lupos stan/DIC_EMEmixSN/DICSN.mixeiv.r")
source("C:/Users/Nelson Lima/Dropbox/Artigo MIXEIV/Simulação e Banco de dados Stan e Jags/Lupos stan/DIC_EMEmixSN/dskewNmix.mult2.r")
source("C:/Users/Nelson Lima/Dropbox/Artigo MIXEIV/Simulação e Banco de dados Stan e Jags/Lupos stan/DIC_EMEmixSN/dSNmixme.r")
source("C:/Users/Nelson Lima/Dropbox/Artigo MIXEIV/Simulação e Banco de dados Stan e Jags/Lupos stan/DIC_EMEmixSN/dSNni2.r")
source("C:/Users/Nelson Lima/Dropbox/Artigo MIXEIV/Simulação e Banco de dados Stan e Jags/Lupos stan/estimate_SN/Stan Model simple SN.r")


fitSN   <- rstan::stan(model_code=stan_codeSN, iter=n.iter,warmup = n.burn, thin = n.thin, chains=1, 
                data = list(N = nrow(Z), xobs = Z[,1], y = Z[,2]),save_dso = FALSE)
nparSN = 11
codaSN = as.data.frame(fitSN)
namesSN = names(codaSN)[1:nparSN]
windows(width=1e+8, height=75e+5)

par(mfrow=c(2,6))
plot.ts(codaSN$beta,main="beta")
plot.ts(codaSN$alpha,main="alpha")
plot.ts(codaSN$omega2_1,main="omega2_1")
plot.ts(codaSN$omega2_0,main="omega2_0")
plot.ts(codaSN$"mu[1]",main="mu[1]")
plot.ts(codaSN$"mu[2]",main="mu[2]")
plot.ts(codaSN$"Delta[1]",main="Delta[1]")
plot.ts(codaSN$"Delta[2]",main="Delta[2]")
plot.ts(codaSN$"gamma2[1]",main="gamma2[1]")
plot.ts(codaSN$"gamma2[2]",main="gamma2[2]")
plot.ts(codaSN$peso,main="peso")

savePlot("Traceplot modelo EMEmixSN",type="jpeg");graphics.off()

dicSN = DICSN.mixeiv(Z,alpha=t(codaSN$"alpha"),beta=t(codaSN$"beta"),mu=t(cbind(codaSN$"mu[1]",codaSN$"mu[2]")),
                     gama=t(cbind(codaSN$"gamma2[1]",codaSN$"gamma2[2]")),Delta =t(cbind(codaSN$"Delta[1]",codaSN$"Delta[2]")),
                     peso=t(cbind(codaSN$"peso",1-codaSN$"peso")),Omega=t(cbind(codaSN$"omega2_0",codaSN$"omega2_1")))

matSN = matrix(c(namesSN,round(cbind(colMeans(codaSN),apply(codaSN,2,"median"),apply(codaSN,2,"sd"),HPDinterval(mcmc(codaSN))),4)[1:nparSN,]),nparSN,6,
               dimnames=list(1:nparSN,c("parameters","mean","median","sd","lower","upper")))

capture.output(as.data.frame(matSN),file = "matriz de estimação- SN.txt")
save(list = ls(),file = "rdata.RData")

# ST case

source("C:/Users/Nelson Lima/Dropbox/Artigo MIXEIV/Simulação e Banco de dados Stan e Jags/Lupos stan/DIC_EMEmixST/DICSTmixeiv.r")
source("C:/Users/Nelson Lima/Dropbox/Artigo MIXEIV/Simulação e Banco de dados Stan e Jags/Lupos stan/DIC_EMEmixST/dskewTmix.mult2.r")
source("C:/Users/Nelson Lima/Dropbox/Artigo MIXEIV/Simulação e Banco de dados Stan e Jags/Lupos stan/DIC_EMEmixST/dSTmixme.r")
source("C:/Users/Nelson Lima/Dropbox/Artigo MIXEIV/Simulação e Banco de dados Stan e Jags/Lupos stan/DIC_EMEmixST/dstni2.r")
source("C:/Users/Nelson Lima/Dropbox/Artigo MIXEIV/Simulação e Banco de dados Stan e Jags/Lupos stan/estimate_ST/Stan Model simple ST.r")

initf <- function(){list(nu=2)}
fitST   <- rstan::stan(model_code=stan_codeST, iter=n.iter,warmup = n.burn, thin = n.thin, chains=1, 
                data = list(N = nrow(Z), xobs = Z[,1], y = Z[,2]),save_dso = FALSE, init=initf)
nparST = 12
codaST = as.data.frame(fitST)
namesST = names(codaST)[1:nparST]
windows(width=1e+8, height=75e+5)

par(mfrow=c(2,6))
plot.ts(codaST$beta,main="beta")
plot.ts(codaST$alpha,main="alpha")
plot.ts(codaST$omega2_1,main="omega2_1")
plot.ts(codaST$omega2_0,main="omega2_0")
plot.ts(codaST$"mu[1]",main="mu[1]")
plot.ts(codaST$"mu[2]",main="mu[2]")
plot.ts(codaST$"Delta[1]",main="Delta[1]")
plot.ts(codaST$"Delta[2]",main="Delta[2]")
plot.ts(codaST$"gamma2[1]",main="gamma2[1]")
plot.ts(codaST$"gamma2[2]",main="gamma2[2]")
plot.ts(codaST$peso,main="peso")
plot.ts(codaST$nu,main="nu")

savePlot("Traceplot modelo EMEmixST",type="jpeg");graphics.off()

dicST = DICST.mixeiv(Z,alpha=t(codaST$"alpha"),beta=t(codaST$"beta"),mu=t(cbind(codaST$"mu[1]",codaST$"mu[2]")),
                     gama=t(cbind(codaST$"gamma2[1]",codaST$"gamma2[2]")),Delta =t(cbind(codaST$"Delta[1]",codaST$"Delta[2]")),
                     peso=t(cbind(codaST$"peso",1-codaST$"peso")),Omega=t(cbind(codaST$"omega2_0",codaST$"omega2_1")),nu=codaST$"nu")

matST = matrix(c(namesST,round(cbind(colMeans(codaST),apply(codaST,2,"median"),apply(codaST,2,"sd"),HPDinterval(mcmc(codaST))),4)[1:nparST,]),nparST,6,
               dimnames=list(1:nparST,c("parameters","mean","median","sd","lower","upper")))

capture.output(as.data.frame(matST),file = "matriz de estimação- ST.txt")
save(list = ls(),file = "rdata.RData")

# SNC case

source("C:/Users/Nelson Lima/Dropbox/Artigo MIXEIV/Simulação e Banco de dados Stan e Jags/Lupos stan/DIC_EMEmixSCN/DICSNCmixeiv.r")
source("C:/Users/Nelson Lima/Dropbox/Artigo MIXEIV/Simulação e Banco de dados Stan e Jags/Lupos stan/DIC_EMEmixSCN/dskewNCmix.mult2.r")
source("C:/Users/Nelson Lima/Dropbox/Artigo MIXEIV/Simulação e Banco de dados Stan e Jags/Lupos stan/DIC_EMEmixSCN/dsnci2.r")
source("C:/Users/Nelson Lima/Dropbox/Artigo MIXEIV/Simulação e Banco de dados Stan e Jags/Lupos stan/DIC_EMEmixSCN/dSTmixme.r")
source("C:/Users/Nelson Lima/Dropbox/Artigo MIXEIV/Simulação e Banco de dados Stan e Jags/Lupos stan/estimate_SCN/Stan Model simple SCN.r")

#initf <- function(){list(beta=0.4,alpha=0.16,omega2_1=0.08,omega2_0=0.005,mu=c(0,8.6),gamma2=c(0.001,0.001),Delta=c(1.1,-0.07),peso=c(0.88,0.11),tau=0.02,rho=0.2)}
fitSCN   <- rstan::stan(model_code=stan_codeSCN, iter=n.iter,warmup = n.burn, thin = n.thin, chains=1, 
                data = list(N = nrow(Z), xobs = Z[,1], y = Z[,2]),save_dso = FALSE)
nparSCN = 13
codaSCN = as.data.frame(fitSCN)
namesSCN = names(codaSCN)[1:nparSCN]
windows(width=1e+8, height=75e+5)

par(mfrow=c(3,5))
plot.ts(codaSCN$beta,main="beta")
plot.ts(codaSCN$alpha,main="alpha")
plot.ts(codaSCN$omega2_1,main="omega2_1")
plot.ts(codaSCN$omega2_0,main="omega2_0")
plot.ts(codaSCN$"mu[1]",main="mu[1]")
plot.ts(codaSCN$"mu[2]",main="mu[2]")
plot.ts(codaSCN$"Delta[1]",main="Delta[1]")
plot.ts(codaSCN$"Delta[2]",main="Delta[2]")
plot.ts(codaSCN$"gamma2[1]",main="gamma2[1]")
plot.ts(codaSCN$"gamma2[2]",main="gamma2[2]")
plot.ts(codaSCN$peso,main="peso")
plot.ts(codaSCN$tau,main="tau")
plot.ts(codaSCN$rho,main="rho")

savePlot("Traceplot modelo EMEmixSCN",type="jpeg");graphics.off()

dicSCN = DICSNC.mixeiv(Z,alpha=t(codaSCN$"alpha"),beta=t(codaSCN$"beta"),mu=t(cbind(codaSCN$"mu[1]",codaSCN$"mu[2]")),
                       gama=t(cbind(codaSCN$"gamma2[1]",codaSCN$"gamma2[2]")),Delta =t(cbind(codaSCN$"Delta[1]",codaSCN$"Delta[2]")),
                       peso=t(cbind(codaSCN$"peso",1-codaSCN$"peso")),Omega=t(cbind(codaSCN$"omega2_0",codaSCN$"omega2_1")),
                       tau=codaSCN$"tau",rho=codaSCN$"rho")

matSCN = matrix(c(namesSCN,round(cbind(colMeans(codaSCN),apply(codaSCN,2,"median"),apply(codaSCN,2,"sd"),HPDinterval(mcmc(codaSCN))),4)[1:nparSCN,]),nparSCN,6,
               dimnames=list(1:nparSCN,c("parameters","mean","median","sd","lower","upper")))

capture.output(as.data.frame(matSCN),file = "matriz de estimação- SCN.txt")
save(list = ls(),file = "rdata.RData")

# SLL case

source("C:/Users/Nelson Lima/Dropbox/Artigo MIXEIV/Simulação e Banco de dados Stan e Jags/Lupos stan/DIC_EMEmixSSL/DICSSLmixeiv.r")
source("C:/Users/Nelson Lima/Dropbox/Artigo MIXEIV/Simulação e Banco de dados Stan e Jags/Lupos stan/DIC_EMEmixSSL/dskewSSL.r")
source("C:/Users/Nelson Lima/Dropbox/Artigo MIXEIV/Simulação e Banco de dados Stan e Jags/Lupos stan/DIC_EMEmixSSL/dskewSSLmix.mult.r")
source("C:/Users/Nelson Lima/Dropbox/Artigo MIXEIV/Simulação e Banco de dados Stan e Jags/Lupos stan/DIC_EMEmixSSL/dSSLmixme.r")
source("C:/Users/Nelson Lima/Dropbox/Artigo MIXEIV/Simulação e Banco de dados Stan e Jags/Lupos stan/estimate_SSL/Stan Model simple SSL.r")

#initf <- function(){list(beta=0.4,alpha=0.16,omega2_1=0.08,omega2_0=0.005,mu=c(0,8.6),gamma2=c(0.001,0.001),Delta=c(1.1,-0.07),peso=c(0.88,0.11),tau=0.02,rho=0.2)}
fitSSL   <- rstan::stan(model_code=stan_codeSSL, iter=n.iter,warmup = n.burn, thin = n.thin, chains=1, 
                 data = list(N = nrow(Z), xobs = Z[,1], y = Z[,2]),save_dso = FALSE)
nparSSL = 12
codaSSL = as.data.frame(fitSSL)
namesSSL = names(codaSSL)[1:nparSSL]
windows(width=1e+8, height=75e+5)

par(mfrow=c(3,5))
plot.ts(codaSSL$beta,main="beta")
plot.ts(codaSSL$alpha,main="alpha")
plot.ts(codaSSL$omega2_1,main="omega2_1")
plot.ts(codaSSL$omega2_0,main="omega2_0")
plot.ts(codaSSL$"mu[1]",main="mu[1]")
plot.ts(codaSSL$"mu[2]",main="mu[2]")
plot.ts(codaSSL$"Delta[1]",main="Delta[1]")
plot.ts(codaSSL$"Delta[2]",main="Delta[2]")
plot.ts(codaSSL$"gamma2[1]",main="gamma2[1]")
plot.ts(codaSSL$"gamma2[2]",main="gamma2[2]")
plot.ts(codaSSL$peso,main="peso")
plot.ts(codaSSL$nu,main="nu")


savePlot("Traceplot modelo EMEmixSSL",type="jpeg");graphics.off()

dicSSL = DICSSL.mixeiv(Z,alpha=t(codaSSL$"alpha"),beta=t(codaSSL$"beta"),mu=t(cbind(codaSSL$"mu[1]",codaSSL$"mu[2]")),
                       gama=t(cbind(codaSSL$"gamma2[1]",codaSSL$"gamma2[2]")),Delta =t(cbind(codaSSL$"Delta[1]",codaSSL$"Delta[2]")),
                       peso=t(cbind(codaSSL$"peso",1-codaSSL$"peso")),Omega=t(cbind(codaSSL$"omega2_0",codaSSL$"omega2_1")),
                       nu=codaSSL$"nu")

matSSL = matrix(c(namesSSL,round(cbind(colMeans(codaSSL),apply(codaSSL,2,"median"),apply(codaSSL,2,"sd"),HPDinterval(mcmc(codaSSL))),4)[1:nparSSL,]),nparSSL,6,
                dimnames=list(1:nparSSL,c("parameters","mean","median","sd","lower","upper")))

capture.output(as.data.frame(matSSL),file = "matriz de estimação- SSL.txt")

Dic = matrix(as.matrix(cbind(dicN,dicT,dicSL,dicCN,dicSN,dicST,dicSSL,dicSCN)),3,8,dimnames=list(c("DIC.obs","tauD","Llog.hat"),c("N","T","SL","NC","SN","ST","SSL","SNC")))
#Dic = matrix(as.matrix(cbind(dicN,dicT,dicSL,dicCN,dicSN,dicST,dicSSL)),3,7,dimnames=list(c("DIC.obs","tauD","Llog.hat"),c("N","T","NC","SL","SN","ST","SSL")))

capture.output(Dic, file="DIC.txt")

#grafico para artigo

postscript(paste(getwd(),"/Estimate-ST.eps",sep=''),family="Times",width=3.2,height=3.2,paper="a4",pointsize=8)
par(mar=c(3.9,4.2,1.4,0.5),mgp=c(2.4, 1, 0))
op <- par(mfrow = c(1,1))
with(codaST, {
  plot.ts(codaST$beta,main=expression(beta),xlab="",ylab="")
  plot.ts(codaST$alpha,main=expression(alpha),xlab="",ylab="")
  plot.ts(codaST$omega2_0,main=expression(omega[0]^2),xlab="",ylab="")
  plot.ts(codaST$omega2_1,main=expression(omega[1]^2),xlab="",ylab="")
  plot.ts(codaST$"mu[1]",main=expression(mu[1]),xlab="",ylab="")
  plot.ts(codaST$"mu[2]",main=expression(mu[2]),xlab="",ylab="")
  plot.ts(codaST$"Delta[1]",main=expression(Delta[1]),xlab="",ylab="")
  plot.ts(codaST$"Delta[2]",main=expression(Delta[2]),xlab="",ylab="")
  plot.ts(codaST$"gamma2[1]",main=expression(gamma[1]^2),xlab="",ylab="")
  plot.ts(codaST$"gamma2[2]",main=expression(gamma[2]^2),xlab="",ylab="")
  plot.ts(codaST$peso,main="p",xlab="",ylab="")
  plot.ts(codaST$nu,main=expression(nu),xlab="",ylab="")
})
mtext(" ",outer=TRUE, font = par("font.main"), cex = 0.5)
par(op)
dev.off()

save(list = ls(),file = "rdata.RData")
