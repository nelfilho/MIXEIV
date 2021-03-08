pd.solve <- function(x, silent=FALSE, log.det=FALSE)
{
  if(is.null(x)) return(NULL)
  if(any(is.na(x)))
  {if(silent) return (NULL) else stop("NA's in x") } 
  if(length(x) == 1) x <- as.matrix(x)
  if(!(inherits(x, "matrix"))) 
  {if(silent) return(NULL) else stop("x is not a matrix")}
  if(max(abs(x - t(x))) > .Machine$double.eps) 
  {if(silent) return (NULL) else stop("x appears to be not symmetric") } 
  x <- (x + t(x))/2
  u <- try(chol(x, pivot = FALSE), silent = silent)
  if(inherits(u, "try-error")) {
    if(silent) return(NULL) else
      stop("x appears to be not positive definite") }
  inv <- chol2inv(u)
  if(log.det) attr(inv, "log.det") <- 2 * sum(log(diag(u)))
  dimnames(inv) <- rev(dimnames(x))
  return(inv)
}

.onLoad <- function(library, pkg)
{ 
  library.dynam("mnormt", pkg, library)
  invisible()
}


dmnorm <- function(x, mean=rep(0,d), varcov, log=FALSE)
{
  d <- if(is.matrix(varcov)) ncol(varcov) else 1
  if(d==1) return(dnorm(x, mean, sqrt(varcov), log=log))
  x <- if (is.vector(x)) t(matrix(x)) else data.matrix(x)
  if(ncol(x) != d) stop("mismatch of dimensions of 'x' and 'varcov'")
  if(is.matrix(mean)) { if ((nrow(x) != nrow(mean)) || (ncol(mean) != d))
    stop("mismatch of dimensions of 'x' and 'mean'") }
  if(is.vector(mean)) mean <- outer(rep(1, nrow(x)), as.vector(matrix(mean,d)))
  X  <- t(x - mean)
  conc <- pd.solve(varcov, log.det=TRUE)
  Q <- colSums((conc %*% X)* X)
  log.det <- attr(conc, "log.det")
  logPDF <- as.vector(Q + d*logb(2*pi) + log.det)/(-2)
  if(log) logPDF else exp(logPDF)
}

dmt <- function (x, mean=rep(0,d), S, df = Inf, log = FALSE) 
{
  if (df == Inf)  return(dmnorm(x, mean, S, log = log))
  d <- if(is.matrix(S)) ncol(S) else 1
  if (d==1) {
    y <- dt((x-mean)/sqrt(S), df=df, log=log)
    if(log) y <- (y - 0.5*logb(S)) else y <- y/sqrt(S)
    return(y)
  }
  x <- if (is.vector(x)) t(matrix(x)) else data.matrix(x)
  if (ncol(x) != d) stop("mismatch of dimensions of 'x' and 'varcov'")
  if (is.matrix(mean)) {if ((nrow(x) != nrow(mean)) || (ncol(mean) != d))
    stop("mismatch of dimensions of 'x' and 'mean'") }
  if(is.vector(mean)) mean <- outer(rep(1, nrow(x)), as.vector(matrix(mean,d)))
  X  <- t(x - mean)
  S.inv <- pd.solve(S, log.det=TRUE)
  Q <- colSums((S.inv %*% X) * X)
  logDet <- attr(S.inv, "log.det")
  logPDF <- (lgamma((df + d)/2) - 0.5 * (d * logb(pi * df) + logDet)
             - lgamma(df/2) - 0.5 * (df + d) * logb(1 + Q/df))
  if(log) logPDF else exp(logPDF)
}
