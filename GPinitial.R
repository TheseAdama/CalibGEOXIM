# Installation package
install.packages("C:/Users/barrya/Desktop/SVDGP.tar.gz",
                 repos = NULL, type = "source")

library(SVDGP)
library(ggplot2)
library(patchwork)
library(tidyr)
library(reshape2)
library(parallel)
library(tictoc)

# Chargement
THETA <- readRDS("C:/Users/barrya/Desktop/CalibGEOXIM/Datas/THETA.rds")
Fsim <- readRDS("C:/Users/barrya/Desktop/CalibGEOXIM/Datas/Fsim.rds")
Dinit <- readRDS("C:/Users/barrya/Desktop/CalibGEOXIM/Datas/Dinit.rds")
Finit <- readRDS("C:/Users/barrya/Desktop/CalibGEOXIM/Datas/Finit.rds")
DD <- readRDS("C:/Users/barrya/Desktop/CalibGEOXIM/Datas/DD.rds")
XX <- readRDS("C:/Users/barrya/Desktop/CalibGEOXIM/Datas/XX.rds")
SimCorresp <- function(D, XX, THETA, Fsim, Nt=461){
  n = nrow(D)
  inds = indt = indx = c()
  for (i in 1:n) {
    ii = which(THETA[, 1]==as.numeric(D[i, 1]) & THETA[, 2]==as.numeric(D[i, 2])
               & THETA[, 3]==as.numeric(D[i, 3])& THETA[, 4]==as.numeric(D[i, 4]))

    jj = which(XX[, 1]==as.numeric(D[i, 5])& XX[, 2]==as.numeric(D[i, 6]))

    indt = c(indt, ii)
    indx = c(indx, jj)
    # cat("Iteration:", i, "\n")
  }
  FD = matrix(0, ncol=Nt, nrow=n)
  for(i in 1:n){
    FD[i, ] <- Fsim[indx[i], , indt[i]]
  }
  return(FD)
}

# Emulateur GP initial
Finit = 1e5*SimCorresp(Dinit, XX, THETA, Fsim, Nt=461)
K=4
model <- svdgppmodel(D=Dinit, FD=Finit ,
                    formula = ~.,
                    K = K, spcovtype ="matern5_2",
                    lower=1e-2,
                    typekrig = "SK", multistart = 20,
                    control=list(maxit=1000))

cat("Erreur relative d'approximation", model$Err)
model$GP[[1]]@covariance
model$GP[[2]]@covariance
model$GP[[3]]@covariance
model$GP[[4]]@covariance

# Coefficient de prédictivité
Npred = 5000
aa=sample(1:nrow(DD), Npred)
Dpred=matrix(DD[aa,], nrow=Npred)
FDpred = 1e5*SimCorresp(Dpred, XX, THETA, Fsim, Nt=461)
YY = svdgppredict(model, Dpred, computecov = FALSE)
Mn = YY$Mn
Q2 = 1 - (sum((FDpred - Mn)^ 2)/sum((FDpred - mean(Finit))^2))
print(100*Q2)

# Coefficient de predictivité pour chaque IndGP
U=(FDpred-colMeans(FDpred))%*%model$V
Q2 =rep(0, K)
for(k in 1:K){
  Q2[k] = 1 - (sum((U[,k] - YY$Coefpred[,k])^ 2)/sum((U[,k] - mean(U[,k]))^2))
}
print(100*Q2)
setwd("C:/Users/barrya/Desktop/CalibGEOXIM/Datas")
save(model, file = "modelm0.RData")

# Amélioration par planification séquentielle
Dadd = Dinit
FD = 1e5*Finit
M = 125
t1 = tic()
for(m in 1:M){
  Dnew = matrix(0, ncol=6, nrow=K)
  for(k in 1:K){
    R <- apply(DD, 1, function(u) {
      r = svdgppredict(model, matrix(u, nrow=1))
      o = as.numeric(r$spcov[[k]])
      return(o)
    })
    istar = which.max(R)
    Dnew[k,] = DD[istar,]
    DD = DD[-istar,]
  }
  FDnew = 1e5*SimCorresp(Dnew, XX, THETA, Fsim, Nt=461)
  model <- updatesvdgp(model, Dnew, FDnew)
  Dadd = rbind(Dadd, Dnew)
  cat("Iteration : ", m, "\n")
}
t2 = toc()
temp = t2$toc -t1

# Q2
Dpred = Dxtr[sample(1:nrow(Dxtr), 1000),]
FDpred = 1e5*SimCorresp(Dpred, XX, THETA, Fsim, Nt=461)
YY = svdgppredict(model=model, Dpred)
Mn = YY$Mn
Qg2 = 1 - (sum((FDpred - Mn)^ 2)/sum((FDpred - mean(FDpred))^2))
print(100*Qg2)

# Save
setwd("C:/Users/barrya/Desktop/CalibGEOXIM/Datas")
saveRDS(temp, file = "tempsvark.rds")
saveRDS(Dadd, file = "Dadd.rds")
saveRDS(DD, file = "DDr.rds")
save(model, file = "modelm1.RData")


