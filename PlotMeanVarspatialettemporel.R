XX <- readRDS("C:/Users/barrya/Desktop/CalibGEOXIM/Datas/XX.rds")
THETA <- readRDS("C:/Users/barrya/Desktop/CalibGEOXIM/Datas/THETA.rds")
Fsim <- readRDS("C:/Users/barrya/Desktop/CalibGEOXIM/Datas/Fsim.rds")
load("C:/Users/barrya/Desktop/CalibGEOXIM/Datas/model0.RData")
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

library(viridis)
library(SVDGP)
Nx = 5000
Nt = 461
Mx = Vx = rep(0, Nx)
for(i in 1:Nx){
  x = XX[i, ]
  Dpred = cbind(THETA, matrix(rep(x, 100), ncol=2, byrow = TRUE))
  FDpred = 1e5*SimCorresp(Dpred, XX, THETA, Fsim, Nt=461)
  YY = svdgppredict(model, Dpred)
  Mx[i] = mean(YY$Mn)
  Vx[i] = sqrt(mean((YY$Mn - mean(YY$Mn))**2))
  cat("Itération : ", i, "\n")
}

Mx = matrix(Mx, nrow = 50, ncol = 100, byrow = FALSE)
Vx = matrix(Vx, nrow = 50, ncol = 100, byrow = FALSE)
x = seq(0, 10, len=50)
y = seq(0, 20, len=100)

setwd("C:/Users/barrya/Desktop/CalibGEOXIM/Figures")
par(mfrow = c(1,1))
png(filename = "Mxhat.png", width = 1800, height = 1900, res = 300)
image(x, y, Mx, col = viridis(100),
      xlab = "X (m)", ylab = "Y (m)",
      main =expression(10^{5}~hat(Mx)), asp = 1)
filled.contour(x, y, Mx, color.palette = viridis,
               plot.title = title(main = expression(10^{5}~hat(Mx)),
                                  xlab = "x", ylab = "y"),
               plot.axes = {axis(1); axis(2)},
               key.title = expression(10^{5}~hat(Mx)))

dev.off()
png(filename = "Vxhat.png", width = 1800, height = 1900, res = 300)
image(x, y, Vx, col = viridis(100), xlab = "X (m)", ylab = "Y (m)",
      main = expression(10^{5}~hat(Vx)), asp = 1)
filled.contour(x, y, Vx, color.palette = viridis,
               plot.title = title(main = expression(10^{5}~hat(Vx)),
                                  xlab = "x", ylab = "y"),
               plot.axes = {axis(1); axis(2)},
               key.title = expression(10^{5}~hat(Vx)))
dev.off()

setwd("C:/Users/barrya/Desktop/CalibGEOXIM/Datas")
saveRDS(Mx, "Mxmodel.rds")
saveRDS(Vx, "Vxmodel.rds")

