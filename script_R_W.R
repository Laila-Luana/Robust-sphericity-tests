library(robustbase)
EstCov <- function(X)
{
  C <- covComed(X)$cov
  S <- var(X)
  return(list(C = C, S = S))
}

WAs <- function(n, Ca1)
{
  p <- nrow(Ca1)
  chi2 <- sum(diag((Ca1 - diag(p)) %*% (Ca1 - diag(p)))) * n / 2
  chi2 <- chi2 - (sum(diag(Ca1)) / p)^2 * p^2 / 2 + p^2 / 2
  p.valor <- 1 - pchisq(chi2, p * (p + 1) / 2)
  return(list(chi2 = chi2, valor.p = p.valor))
}  

WMC <-  function(N = 1000, n, Ca)
{
  chi2C <- WAs(n, Ca$C)$chi2
  chi2S <- WAs(n, Ca$S)$chi2
  library(MASS)
  p <- nrow(Ca$C)
  mu <- rep(0, times  = p)
  Sigma <- diag(p)
  for (i in 1:N)
  {
    X <- mvrnorm(n, mu, Sigma)
    Ca <- EstCov(X)
    chi2C <- c(chi2C, WAs(n, Ca$C)$chi2)
    chi2S <- c(chi2S, WAs(n, Ca$S)$chi2)
  }
  valorC.p <- sum(chi2C >= chi2C[1]) / (N + 1)
  valorS.p <- sum(chi2S >= chi2S[1]) / (N + 1)
  return(list(Chi2C = chi2C[1], valorC.p = valorC.p,
              Chi2S = chi2S[1], valorS.p = valorS.p))
}  

library(MASS)
rNCM <- function(n, delta , mu, mu2,  Sigma,  Sigma2)
{
  u <- runif(n)
  p <- nrow(Sigma)
  n1 <- length(u[u <= delta])
  if (n1 < 1) n1 <- 1
  n2 <- n - n1
  X <- matrix(0, n, p)
  X[u <= delta, ] <- mvrnorm(n1, mu, Sigma)
  if (n2 > 0) X[u > delta, ] <- mvrnorm(n2, mu2, Sigma2)
  return(X)
}

library(mvtnorm)
simMCpgen <- function(op, opH0, n, p, rho, delta, N, NMC, k)
{
  rej <- matrix(0, 4, 3)
  rownames(rej) <- c("WAs", "WAsR","WMC", "WMCR")
  colnames(rej) <- c("0,10","0,05","0,01")
  N1R <- 1.0 / N
  mu <- rep(0, times  = p)
  if (opH0 == 0) Sigma <- diag(p) else
    if (opH0 == 1) Sigma <- 2 * diag(p) else
      if (opH0 == 2) Sigma <- diag(runif(p,1, 2)) else
        if (opH0 == 3) Sigma <- (1 - rho) * diag(p) + rho * matrix(1, p, p)
  if (op == 2)
  {
    Delta <- (k - delta) / (1 - delta)
    Sigma2 <- Delta * Sigma
    mu2 <- rep(0, times = p)
  }
  for (i in 1:N)
  {    
    if (op == 1) X <- mvrnorm(n, mu, Sigma) else
      if (op == 2) X <- rNCM(n, delta, mu, mu2, Sigma, Sigma2)
      Ca <- EstCov(X)
      C1 <- Ca$C
      S1 <- Ca$S
      res3 <- WAs(n, S1)
      res4 <- WAs(n, C1)
      res6 <- WMC(NMC, n, Ca)
      if (res3$valor.p <= 0.10)  rej[1,1] <- rej[1,1] + N1R
      if (res3$valor.p <= 0.05)  rej[1,2] <- rej[1,2] + N1R
      if (res3$valor.p <= 0.01)  rej[1,3] <- rej[1,3] + N1R
      if (res4$valor.p <= 0.10)  rej[2,1] <- rej[2,1] + N1R
      if (res4$valor.p <= 0.05)  rej[2,2] <- rej[2,2] + N1R
      if (res4$valor.p <= 0.01)  rej[2,3] <- rej[2,3] + N1R
      if (res6$valorS.p <= 0.10) rej[3,1] <- rej[3,1] + N1R
      if (res6$valorS.p <= 0.05) rej[3,2] <- rej[3,2] + N1R
      if (res6$valorS.p <= 0.01) rej[3,3] <- rej[3,3] + N1R
      if (res6$valorC.p <= 0.10) rej[4,1] <- rej[4,1] + N1R
      if (res6$valorC.p <= 0.05) rej[4,2] <- rej[4,2] + N1R
      if (res6$valorC.p <= 0.01) rej[4,3] <- rej[4,3] + N1R
  }
  return(rej)
}  

op <- 2
opH0 <- 0
n <- 4
p <- 8
k <- 1.005
rho <- 0.9
delta <- 0.9
N <- 999
NMC <- 1000

res <- simMCpgen(op, opH0, n, p, rho, delta, N, NMC, k)
res
