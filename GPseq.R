# Chargement package
library(SVDGP)
library(ggplot2)

# Chargement
THETA <- readRDS("C:/Users/barrya/Desktop/CalibGEOXIM/Datas/THETA.rds")
Fsim <- readRDS("C:/Users/barrya/Desktop/CalibGEOXIM/Datas/Fsim.rds")
Dinit <- readRDS("C:/Users/barrya/Desktop/CalibGEOXIM/Datas/Dinit.rds")
Finit <- readRDS("C:/Users/barrya/Desktop/CalibGEOXIM/Datas/Finit.rds")
DD <- readRDS("C:/Users/barrya/Desktop/CalibGEOXIM/Datas/DD.rds")
XX <- readRDS("C:/Users/barrya/Desktop/CalibGEOXIM/Datas/XX.rds")
Xobs <- readRDS("C:/Users/barrya/Desktop/Application/Xopt.rds")
Yobs <- readRDS("C:/Users/barrya/Desktop/Application/Yopt.rds")
load("C:/Users/barrya/Desktop/CalibGEOXIM/Datas/modelm0.RData")

# Definition des fonctions
Unif <- function(K, D){
  R = matrix(0, nrow=K, ncol=nrow(D))
  for(k in 1:nrow(D)){
    R[,k] <- runif(K, min=D[k,1], max = D[k,2])}
  return(R)
}
minSSk = function(Yopt, Xopt, model, Lmc, DomT){
  TT = Unif(Lmc, DomT)
  R  = apply(TT, 1, function(theta){
    D = cbind(matrix(rep(theta, 8), ncol=4, byrow = TRUE), Xopt)
    D = as.matrix(D)
    Ysim =simsvdgp(model, D, Lmc)
    SS = apply(array(c(1:Lmc)), 1, function(l){
      R = norm(1e5*Yopt-Ysim[l,,])**2
      return(R)
    })
    o = mean(SS,na.rm = TRUE)
    return(o)
  })

  o = min(R, na.rm = TRUE)
  return(o)
}
EIk = function(theta, model, Yopt, Xopt, mSSk, Lmc){

  D = cbind(matrix(rep(theta, 8), ncol=4, byrow = TRUE), Xopt)
  D = as.matrix(D)
  Ysim =simsvdgp(model, D, Lmc)
  SS = apply(array(c(1:Lmc)), 1, function(l){
    R = norm(1e5*Yopt-Ysim[l,,])**2
    return(R)
  })
  o = mean(max(mSSk-SS, 0))
  return(o)
}

Svarvk <- function(x, theta, model){
  D = cbind(matrix(theta, nrow=1), matrix(x, nrow=1))
  Yp = svdgppredict(model, Dpred =  matrix(as.numeric(D), nrow=1))
  o = sum(c(Yp$spcov[[1]],Yp$spcov[[2]],Yp$spcov[[3]],Yp$spcov[[4]]))
  return(o)
}

# Planification sequentielle
DomT = matrix(c(0.2, 0.4, 0.2,-12.31,1.6,4.4, 0.45,-10.39), ncol=2)
D <- readRDS("C:/Users/barrya/Desktop/Application/Dadd.rds")
FD = 1e5*SimCorresp(D, XX, THETA, Fsim, Nt=461)
Lmc = 100
M = 500
m = 1
ivec = jvec = c(0)
while(m<=M+1){
  mSSk = minSSk(Yopt, Xopt, model, Lmc, DomT)
  EIv = apply(THETA, 1, function(theta){
    R = EIk(theta, model, Yopt, Xopt, mSSk, Lmc)
    return(R)
  })
  istar = which.max(EIv)
  thetastar = THETA[istar,]
  Svarv = apply(XX, 1, function(x){
    R = Svarvk(x, thetastar, model)
    return(R)
  })
  jstar = which.max(Svarv)
  xstar = XX[jstar,]
  if(any(apply(array(cbind(ivec, jvec)), 1, function(col){all(col == c(istar,jstar))}))){
  }else{
    Dnew = cbind(matrix(thetastar, nrow=1), matrix(xstar, nrow=1))
    Dnew = matrix(as.numeric(Dnew), nrow=1)

    # Mise a jour du modele
    FDnew = 1e5*SimCorresp(Dnew, XX, THETA, Fsim, Nt=461)
    model <- updatesvdgp(model, Dnew, FDnew)
    D = rbind(D, Dnew)
    FD = rbind(FD, FDnew)

    # Mise à jour
    ivec = c(ivec, istar)
    jvec = c(jvec, jstar)
    m = m+1
  }
  cat("Iteration : %d", m)
}

# Save
setwd("C:/Users/barrya/Desktop/CalibGEOXIM/Datas")
saveRDS(D, file = "Dseq.rds")
saveRDS(FD, file = "FDseq.rds")
save(model, file = "modelm2.RData")



