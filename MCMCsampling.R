library(SVDGP)
library(MASS)
# Chargement
load("C:/Users/barrya/Desktop/Application/model2.RData")
THETA <- readRDS("C:/Users/barrya/Desktop/Application/THETA.rds")
Xobs <- readRDS("C:/Users/barrya/Desktop/Application/Xopt.rds")
Yobs <- readRDS("C:/Users/barrya/Desktop/Application/Yopt.rds")

# Tirage uniforme dans un domaine D
Unif <- function(K, D){
  R = matrix(0, nrow=K, ncol=nrow(D))
  for(k in 1:nrow(D)){
    R[,k] <- runif(K, min=D[k,1], max = D[k,2])}
  return(R)
}
# Echantillonnage MCMC
MetropInGibbs <- function(logfd, maxiter, xinit, sig=0.5, Dx=NULL) {
  d <- length(xinit)
  sig = rep(sig, d)
  if(!is.null(Dx)){sig <- apply(Dx, 1, function(u){return(sig*(u[2]-u[1]))})}
  samples <- matrix(NA, nrow = maxiter, ncol = d)
  samples[1, ] <- xinit
  for (i in 2:maxiter) {
    xi <- samples[i - 1, ]
    for (j in 1:d) {
      xprop <- xi
      xprop[j] <- rnorm(1, mean = xi[j], sd = sig[j])
      rr <- logfd(xprop) - logfd(xi)
      if (log(runif(1)) < rr) {
        xi[j] <- xprop[j]
      }

    }
    samples[i, ] <- xi
    cat("ITERATION :", i,"\n")
    if(i%%1000==0){cat("ITERATION :", i,"\n")}
  }

  # samples <- samples[(round(0.2*maxiter) + 1):maxiter, ]
  # samples <- samples[seq(1, nrow(samples), by=3), ]
  return(samples)
}

# Densité a posteriori
Posteriordf <- function(theta, Xobs, Yobs, sigeps,
                        model){
  Yobs = matrix(Yobs, nrow=8)
  Dp = cbind(matrix(rep(theta, nrow(Xobs)), ncol=4, byrow = TRUE), Xobs)
  Yp = svdgppredict(model, Dpred=Dp,computecov = TRUE)
  Mm = 1e-5*Yp$Mn
  Km = Yp$Kn
  Km = 1e-10*Km + diag(sigeps, nrow(Km))
  V = matrix(as.vector(Yobs - Mm), nrow=1)
  SSm = as.numeric(V%*%solve(Km)%*%t(V))
  P = -0.5*SSm - 0.5*log(det(Km)+1e-10)
  return(P)
}

# Log de la dap
logfd = function(u){
  if((sum(u<D[,2])==5)&(sum(u>D[,1])==5)){
    o = Posteriordf(theta=u[1:4], Xobs, Yobs, sigeps=u[5], model)
  }else{o=-Inf}
  return(o)
}
# Domaine des parametres
D = matrix(c(0.2, 0.4, 0.2, -12.31, 2.07e-14,
             1.6, 4.4, 0.45,-10.39, 3.32e-13), ncol=2)

# Lancement MCMCs
setwd("C:/Users/barrya/Desktop/CalibGEOXIM/Datas")
for(i in 1:10){
  xinit <- Unif(1, D)
  R <- MetropInGibbs(logfd = logfd, maxiter = 100, xinit = xinit,
                     sig = 0.175, Dx = D)
  nom <- paste0("MCMC_chain_", i, ".rds")
  saveRDS(R, file = nom)
}
