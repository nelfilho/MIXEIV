##FUNÇÕES NECESSÁRIAS -----------------------------------------------------------------------------------------------------------------------------
##-------------------------------------------------------------------------------------------------------------------------------------------------
matrix.sqrt <- function(A) {
  sva <- svd(A)
  if (min(sva$d)>=0)
    Asqrt <- t(sva$v %*% (t(sva$u) * sqrt(sva$d)))
  else
    stop("Matrix square root is not defined")
  return(Asqrt)
}

dmvt.ls <- function(y, xi, Gama, Delta, nu){
# Calcula a densidade do vetor x para o modelo ST multivariado sob a reparametrização de Cristian
# y: deve ser uma matrix onde cada linha tem um vetor de dados multivariados de dimensão ncol(y) = p. nrow(y) = tamanho da amostra, se caso y seja vetor deve-se colocar em forma de matrix linha
# mu, Delta: devem ser do tipo vetor de mesma dimensão igual a ncol(y) = p
# Gama: Matrix p x p
  
  n <- nrow(y)
  p <- ncol(y)

  Omega  <- Gama + Delta%*%t(Delta)
  lambda <- (solve(matrix.sqrt(Omega))%*%Delta) / as.numeric((1-t(Delta)%*%solve(Omega)%*%Delta))^(1/2)

  denst <- (gamma((p+nu)/2)/(gamma(nu/2)*pi^(p/2)))*nu^(-p/2)*det(Omega)^(-1/2)*(1 + mahalanobis(y, xi, Omega)/nu)^(-(p+nu)/2)
  dens  <- 2*(denst)*pt(sqrt((p + nu)/(mahalanobis(y, xi, Omega) + nu))*apply( matrix(rep(t(lambda)%*%solve(matrix.sqrt(Omega)),n), n, p, byrow = TRUE)*(y - matrix(rep(xi, n), n, p, byrow = TRUE)), 1, sum  )    , df = nu + p   )
  return(dens)
}

bdiag <- function(x){
# Essa função constroi uma matriz bloco diagonal através de uma lista de matrizes
# x é uma lista de matrizes
     if(!is.list(x)) stop("x not a list")
     n <- length(x)
     if(n==0) return(NULL)
     x <- lapply(x, function(y) if(length(y)) as.matrix(y) else
stop("Zero-length component in x"))
     d <- array(unlist(lapply(x, dim)), c(2, n))
     rr <- d[1,]
     cc <- d[2,]
     rsum <- sum(rr)
     csum <- sum(cc)
     out <- array(0, c(rsum, csum))
     ind <- array(0, c(4, n))
     rcum <- cumsum(rr)
     ccum <- cumsum(cc)
     ind[1,-1] <- rcum[-n]
     ind[2,] <- rcum
     ind[3,-1] <- ccum[-n]
     ind[4,] <- ccum
     imat <- array(1:(rsum * csum), c(rsum, csum))
     iuse <- apply(ind, 2, function(y, imat) imat[(y[1]+1):y[2],
(y[3]+1):y[4]], imat=imat)
     iuse <- as.vector(unlist(iuse))
     out[iuse] <- unlist(x)
     return(out)
}

prodarray = function(a,b){
# Essa função faz a multiplicação de cada termo do vetor a pela sua correspondente matriz no array b
# a é um vetor e b é um array de matrizes
c = array(,dim(b))
for(i in 1:length(a)){
c[,,i] = a[i]*b[,,i]
}
return(c)
}

colVars <- function(x, na.rm=FALSE, dims=1, unbiased=TRUE, SumSquares=FALSE,twopass=FALSE) {
#Semelhante ao colMeans e colSums só que calcula variâncias
  if (SumSquares) return(colSums(x^2, na.rm, dims))
  N <- colSums(!is.na(x), FALSE, dims)
  Nm1 <- if (unbiased) N-1 else N
  if (twopass) {x <- if (dims==length(dim(x))) x - mean(x, na.rm=na.rm) else
                     sweep(x, (dims+1):length(dim(x)), colMeans(x,na.rm,dims))}
  (colSums(x^2, na.rm, dims) - colSums(x, na.rm, dims)^2/N) / Nm1
}

d.mixedST <- function (x, pi1, mu, gamma2, Delta, nu){
# Calcula a densidade mistura de ST (univariada) da variável x associada aos vetores de parâmetros
# mu, Delta, gamma2 entrando como vetores de parâmetros reparametrizados (reparametrização de Cristian)
    g <- length(pi1)
    dens <- 0
    for (j in 1:g) dens <- dens + pi1[j] * dt.ls(x, mu[j], gamma2[j],
        Delta[j], nu)
    return(dens)
}

dt.ls <- function (x, loc = 0, gamma2 = 1, Delta = 1, nu = 4){
# Calcula a densidade ST (univariada) da variável x associada aos parâmetros
# mu, Delta, gamma2 entrando como parâmetros reparametrizados (reparametrização de Cristian)
    shape  <- Delta/sqrt(gamma2)
    sigma2 <- gamma2 + Delta^2
    d <- (x - loc)/sqrt(sigma2)
    dens <- 2 * dt(d, df = nu) * pt(sqrt((1 + nu)/(d^2 + nu)) * d * shape, 1 + nu)/sqrt(sigma2)
    return(dens)
}

d.mixedSN <- function (x, pi1, mu, gamma2, Delta){
# Calcula a densidade mistura de SN (univariada) da variável x associada aos vetores de parâmetros
# mu, Delta, gamma2 entrando como vetores de parâmetros reparametrizados (reparametrização de Cristian)
    g <- length(pi1)
    dens <- 0
    for (j in 1:g) dens <- dens + pi1[j] * dSN(x, mu[j], gamma2[j],
        Delta[j])
    return(dens)
}

dSN <- function (y, mu = 0, gamma2 = 1, Delta = 1){
# Calcula a densidade SN (univariada) da variável y associada aos parâmetros
# mu, Delta, gamma2 entrando como parâmetros reparametrizados (reparametrização de Cristian)
    shape  <- Delta/sqrt(gamma2)
    sigma2 <- gamma2 + Delta^2
    dens <- 2 * dnorm(y, mu, sqrt(sigma2)) * pnorm(shape * ((y -
        mu)/sqrt(sigma2)))
    return(dens)
}

dmvSN <- function (y, mu,  Gama, Delta){
# Calcula a densidade do vetor x para o modelo SN multivariado sob a reparametrização de Cristian
# y: deve ser uma matrix onde cada linha tem um vetor de dados multivariados de dimensão ncol(y) = p. nrow(y) = tamanho da amostra, se caso y seja vetor deve-se colocar em forma de matrix linha
# mu, Delta: devem ser do tipo vetor de mesma dimensão igual a ncol(y) = p
# Gama: Matrix p x p
    n <- nrow(y)
    p <- ncol(y)

    Omega  <- Gama + Delta%*%t(Delta)
    lambda <- (solve(root.matrix(Omega))%*%Delta) / as.numeric((1-t(Delta)%*%solve(Omega)%*%Delta))^(1/2)

    #Omega  <- Gama + tcrossprod(Delta,Delta)
    #lambda <- matrix(solve(matrix.sqrt(Omega),Delta)/c(sqrt(1-crossprod(Delta,solve(Omega,Delta)))),nrow=p,ncol=1)

    dens <- 2 * dmnorm(y, mu, Omega) * pnorm(apply(matrix(rep(t(lambda) %*%
        solve(root.matrix(Omega)), n), n, p, byrow = TRUE) *
        (y - matrix(rep(mu, n), n, p, byrow = TRUE)), 1, sum))
    return(dens)
}




pred=function(x,peso,mu,gamma2,Delta,nu,dens="ST"){
# Calcula a fução preditiva para o modelo IEVmixST e IEVmixSN univariado
# x é o vetor de pontos cuja densidade será calculada
# Todos os vetores tem que entrar em forma de matriz linha onde encada linha tomos o l-ésimo vetor onde l=1,...,namost (tamanho da amostra mcmc)
  G = ncol(mu)
  med=matrix(,nrow(mu),length(x))
  for(j in 1:length(x)){
    for(i in 1:nrow(mu)){
      if(dens=="ST"){
        med[i,j]= d.mixedST(x[j],peso[,i],mu[,i],gamma2[,i],Delta[,i],nu[i])
      }
      if(dens=="SN"){
        med[i,j]= d.mixedSN(x[j],peso[,i],mu[,i],gamma2[,i],Delta[,i])
      }
      if(dens=="N"){
        med[i,j]= d.mixedSN(x[j],peso[,i],mu[,i],gamma2[,i],rep(0,G))
      }
      if(dens=="T"){
        med[i,j]= d.mixedST(x[j],peso[,i],mu[,i],gamma2[,i],rep(0,G),nu[i])
      }
      if(dens=="SS"){
        med[i,j]= d.mixedSS(x[j],peso[,i],mu[,i],gamma2[,i],Delta[,i],nu[i])
      }
    }
  }
  return(colMeans(med))
}


d.mixedmvSN <- function (y, pi1, mu, Sigma, lambda) {
# Calcula a densidade mistura de ST segundo a reparametrização de Cristian
# y: deve ser uma matrix onde cada linha tem um vetor de dados multivariados de dimensão ncol(y) = p. nrow(y) = tamanho da amostra, se caso y seja vetor deve-se colocar em forma de matrix linha
# mu, Lambda: parâmetros de locação e forma, devem ser do tipo lista onde o j-ésimo elemento da lista corresponde aos parâmetros da população j da mistura
# Omega: parâmetro de escala, deve ser do tipo lista onde o j-ésimo elemento da lista corresponde aos parâmetros da população j da mistura onde Omega é uma Matrix p x p
# pi1: são os pesos e devem ser do tipo vetor
    g <- length(pi1)
    dens <- 0
    for (j in 1:g) dens <- dens + pi1[j] * dmvSN(y, mu[[j]], 
        Sigma[[j]], lambda[[j]])
    return(dens)
}


d.mixedmvST <- function (y, pi1, mu, Sigma, lambda, nu){
# Calcula a densidade mistura de ST segundo a reparametrização de Cristian
# y: deve ser uma matrix onde cada linha tem um vetor de dados multivariados de dimensão ncol(y) = p. nrow(y) = tamanho da amostra, se caso y seja vetor deve-se colocar em forma de matrix linha
# mu, Lambda: parâmetros de locação e forma, devem ser do tipo lista onde o j-ésimo elemento da lista corresponde aos parâmetros da população j da mistura
# Omega: parâmetro de escala, deve ser do tipo lista onde o j-ésimo elemento da lista corresponde aos parâmetros da população j da mistura onde Omega é uma Matrix p x p
# pi1: são os pesos e devem ser do tipo vetor
    g <- length(pi1)
    dens <- 0
    for (j in 1:g) dens <- dens + pi1[j] * dmvt.ls(y, mu[[j]], 
        Sigma[[j]], lambda[[j]], nu)
    return(dens)
}


DIC =function(Z,peso,mu,gamma2,Delta,alpha,beta,Omega,nu,Dens="ST"){
# Calcula o critério de seleção de modelos DIC para os modelos IEVmixST (Dens=ST) e IEVmixSN (Dens=SN)
# Z é a matriz de dados com a variável regressora na primeira coluna e a matrix resposta nas demais
# Todos os vetores tem que entrar em forma de matriz linha onde encada linha tomos o l-ésimo vetor onde l=1,...,namost (tamanho da amostra mcmc)

med1 = med2 = vector()
p = ncol(Z)
G = ncol(mu)
for(i in 1:nrow(Z)){
dens =  logdens = vector()
for(l in 1:nrow(mu)){

xi = Sigma = Lambda = list()

for(j in 1:G){
xi[[j]]     = c(0,alpha[l,]) + mu[l,j]*c(1,beta[l,])
Sigma[[j]]  = gamma2[l,j]*tcrossprod(c(1,beta[l,])) + diag(Omega[l,])
  if(Dens == "ST" || Dens == "SN"){
    Lambda[[j]] = Delta[l,j]*c(0,beta[l,])
  }else{
    Lambda[[j]] = rep(0,p)
  }
}

if(Dens == "ST" || Dens == "T"){
  dens[l] <- d.mixedmvST(matrix(Z[i,],1,p),peso[l,],xi,Sigma,Lambda,nu[l])
}
if(Dens == "SN" || Dens == "N"){
  dens[l] <- d.mixedmvSN(matrix(Z[i,],1,p),peso[l,],xi,Sigma,Lambda)
}

logdens[l] <- log(dens[l])
}
med1[i] <- mean(logdens)
med2[i] <- mean(dens)
print(i)
}

DIC <- -4*sum(med1)+2*sum(log(med2))    # Deviance Information Criterion
return(DIC)
}


dnormix.mul <- function(y,p,xi,Omega,log=FALSE){

library(mvtnorm)
library(MCMCpack)
k <- length(p)
n <- nrow(y)
m <- ncol(y)
aa <- matrix(0,n,k)
p  <- matrix(p,nrow=k)

  for(j in 1:k){
    aa[,j] <- dmnorm(y,xi[,j],Omega[,,j],log=log)
  }
  f <- aa%*%p
}


dtmix.mul <- function(y,p,xi,Omega,nu){

library(mvtnorm)
library(MCMCpack)
k <- length(p)
n <- nrow(y)
m <- ncol(y)
aa <- matrix(0,n,k)
p  <- matrix(p,nrow=k)

  for(j in 1:k){
    aa[,j] <- dmt(y,xi[,j],Omega[,,j],nu)
  }
  f <- aa%*%p
}

dSN <- function (y, mu = 0, gamma2 = 1, Delta = 1){
# Calcula a densidade SN (univariada) da variável y associada aos parâmetros
# mu, Delta, gamma2 entrando como parâmetros reparametrizados (reparametrização de Cristian)
    shape  <- Delta/sqrt(gamma2)
    sigma2 <- gamma2 + Delta^2
    dens <- 2 * dnorm(y, mu, sqrt(sigma2)) * pnorm(shape * ((y -
        mu)/sqrt(sigma2)))
    return(dens)
}

### Densidade da SSL  ######
dSS2 <- function(y, mu , gamma2, Delta, nu){
# Calcula a densidade SSL (univariada) da variável y associada aos parâmetros
# mu, Delta, gamma2 entrando como parâmetros reparametrizados (reparametrização de Cristian)

  shape  <- Delta/sqrt(gamma2)
  sigma2 <- gamma2 + Delta^2
  
  resp <- vector(mode = "numeric", length = length(y))
  for (i in 1:length(y)) {
    f <- function(u) 2*nu*u^(nu - 1)*dnorm(y[i],mu,sqrt(sigma2/u))*pnorm(u^(1/2)*shape*(sigma2^(-1/2))*(y[i]-mu))
    resp[i] <- integrate(f,0,1)$value
  }
  return(resp)
}

###########    Densidade da Mistura de SSL's   ##################
d.mixedSS <- function(x, pi1, mu, gamma2, Delta, nu){
  # x: é o vetor de dados
  # outros parametros devem ser do tipo vetor c() de dimensão g (qtd de misturas)
  # com parâmetros reparametrizados (reparametrização de Cristian)da função dSS2
  g <- length(pi1)
  dens <- 0
  for (j in 1:g) dens <- dens + pi1[j]*dSS2(x, mu[j], gamma2[j], Delta[j], nu)
  return(dens)
}


