# Função para gerar observações do modelo de regressão linear com com erros nas
# variáveis
# onde a variável latente é uma mistura de normais inversas Gaussianas  
#========================================================================
# Input 
#========================================================================
# n: tamanho amostral
# alpha, beta: vetores de regressão de  dimensão r; 
# Omega: matriz p x p associada à escala  dos erros de observação, onde p=r+1; 
# mu: vetor de  parâmetros associados à locação da variável latente para cada
# sub-população. dim(mu)=G, onde
# G é o número de sub-populações e mu[i] é o parâmetro de locação associado à
# sub-população i;
# lambda: analogamente, é um vetor de dimensão G contendo os parâmetros reguladores
# de assimetria para cada sub-população;
# gama, delta: constantes positivas, parâmetros do fator de escala - que tem
# distribuição Gaussiana inversa
# pii: vetor de dimensão G com os pesos na mistura
#---------------------------------------------------------------------------
#  Início da função
#---------------------------------------------------------------------------
reiv.NIGmix <- function(n,alpha,beta,Omega,mu,lambda,gama,delta,pii){
library(mnormt)
library(SuppDists)
p <- length(alpha) + 1
G <- length(pii)        # número de componentes na mistura 
#-------------------------------------------------------------------------------
x <- matrix(0,n,1)            # vetor que conterá os valores gerados da variável
                              # latente
Z <- matrix(0,n,p)            # vetor que conterá as observações geradas
#===============================================================================
a <- matrix(c(0,alpha))
B1 <- matrix(c(1,beta))
u <- vector()
w <- matrix(sample(G,size=n,replace=TRUE,prob=pii)) # vetor de alocações. Se z[i]=j
                                                 #o indivíduo i é alocado à classe j
nu <- delta/gama
psi <- delta^2 
for (i in 1:n){
u[i] <- rinvGauss(1,nu,psi)
x[i] <- rnorm(1,mu[w[i]]+u[i]*lambda[w[i]],sqrt(u[i]))
Z[i,] <- matrix(rmnorm(1,a+B1*x[i],u[i]*Omega))
}
return(list(Z=Z,w=w,x=x))
    }
    

# Exemplo
mu <- c(-3,1)
#lambda <- c(-0.2,0.1)
#lambda<-c(-0.5,0.4)
lambda<-c(-2,1)
delta <- 0.7
gama <- 1 
pii <- c(0.4,0.6)
Omega <- diag(c(1,1,1)) 
alpha <- c(0.4,0.1)
beta <- c(0.8,0.9)
n <- 200

library(robustbase)    
y.teste <- reiv.NIGmix(n,alpha,beta,Omega,mu,lambda,gama,delta,pii)
#par(mfrow=c(1,5))
#w<-y.teste$w
#plot(density(y.teste$x[w==1]),main='Group 1')
#plot(density(y.teste$x[w==2]),main='Group 2')
#plot(density(y.teste$x),main='latent variable')
#adjbox(y.teste$x[w==1],main='Group 1')
#adjbox(y.teste$x[w==2],main='Group 2')

#Z=y.teste$Z