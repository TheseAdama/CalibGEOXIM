THETA <- readRDS("C:/Users/barrya/Desktop/CalibGEOXIM/Datas/THETA.rds")
Fsim <- readRDS("C:/Users/barrya/Desktop/CalibGEOXIM/Datas/Fsim.rds")
XX <- readRDS("C:/Users/barrya/Desktop/CalibGEOXIM/Datas/XX.rds")

# Package
library(parallel)
n_cores <- detectCores() - 1

# Maillage (theta,x)
Nx = 5000
Dxt = matrix(0, nrow=500000, ncol=6)
for(i in 1:100){
  Dxt[((i-1)*Nx+1):(Nx*i),1:4] = matrix(rep(as.numeric(THETA[i, ]), each = Nx), ncol = 4)
  Dxt[((i-1)*Nx+1):(Nx*i),5:6] =XX
}
saveRDS(Dxt, file = "Dxt.rds")
# Distance min
dMm <- function(Mat) {
  dm <- as.matrix(dist(Mat))
  diag(dm) <- Inf
  o <- min(dm)
  return(o)
}

# Plan Maximin
M = 500
DD  = Dxt
Dinit = matrix(0, ncol=6, nrow=M)
i1 = sample(1:500000,1)
Dinit[1,] = as.numeric(DD[i1,])
DD = DD[-i1,]
for(i in 2:M){
  dMmvec = apply(DD, 1, function(u){
    Mat = rbind(Dinit[c(1:(i-1)),],u)
    return(dMm(Mat))
  })

  istar = which.max(dMmvec)
  Dinit[i,]  = as.numeric(DD[istar, ])
  DD = DD[-istar, ]
  cat("iteration = ", i, "\n")
}
dim(Dinit)

# Simulation correspondante à un design D
SimCorresp <- function(D, XX, THETA, Fsim, Nt=461){
  n = nrow(D)
  inds = indtheta = indx = c()
  for (i in 1:n) {
    ii = which(THETA[, 1]==as.numeric(D[i, 1]) & THETA[, 2]==as.numeric(D[i, 2])
               & THETA[, 3]==as.numeric(D[i, 3])& THETA[, 4]==as.numeric(D[i, 4]))

    jj = which(XX[, 1]==as.numeric(D[i, 5])& XX[, 2]==as.numeric(D[i, 6]))

    indtheta = c(indtheta, ii)
    indx = c(indx, jj)
    cat("Iteration:", i, "\n")
  }
  FD = matrix(0, ncol=Nt, nrow=n)
  for(i in 1:n){
    FD[i, ] <- Fsim[indx[i], , indtheta[i]]
  }
  return(FD)
}

Finit = SimCorresp(Dinit, XX, THETA, Fsim, Nt=461)

# Sauvegarde
saveRDS(Dinit, file = "Dinit.rds")
saveRDS(Finit, file = "Finit.rds")
saveRDS(DD, file = "DD.rds")
saveRDS(XX, file = "XX.rds")

