library(ggplot2)
library(SVDGP)

XX <- readRDS("C:/Users/barrya/Desktop/Application/XX.rds")
load("C:/Users/barrya/Desktop/CalibGEOXIM/Datas/model0.RData")

# Plan uniforme sur un domaine D
Unif <- function(K, D){
  R = matrix(0, nrow=K, ncol=nrow(D))
  for(k in 1:nrow(D)){
    R[,k] <- runif(K, min=D[k,1], max = D[k,2])
  }
  return(R)
}

# Distance minimale
dmin <- function(Mat) {
  dm <- as.matrix(dist(Mat))
  diag(dm) <- Inf
  return(min(dm))
}
# Variation approchée du code

CVmodel <- function(X, model, DomT, K=1000){
  set.seed(123)
  Theta = Unif(K, D=DomT)
  R = rep(0, length=nrow(X))
  for(i in 1:nrow(X)){
    Dpred = cbind(Theta, matrix(rep(as.numeric(X[i,]),K),
                                nrow=K, byrow = TRUE))
    Yp = svdgppredict(model, Dpred)
    R[i] = mean((Yp$Mn - mean(Yp$Mn))**2)
  }
  o = sum(R)
  return(o)
}

# RS epsvec1, epsvec2
tirage <- function(X,lambda = 1e-2){
  xij = matrix(0, nrow = 15, ncol = 2)
  pij = rep(0, 15)
  c = 1
  for(i in 2:6) {
    for(j in (i+1):7) {
      xij[c, ] =  c(i, j)
      dij =  norm(X[i,] - X[j,], "2")
      pij[c] =  1 / (dij + lambda)
      c = c + 1
    }
  }
  pij = pij / sum(pij)
  ii = sample(1:15, 1, prob = pij)
  xx = xij[ii,]
  jj = sample(1:2, 1, prob = c(0.5,0.5))
  ind = xx[jj]
  return(ind)
}

RSoptim <- function(f, Xinit, K=1000, alpha=0.99, T0 = 1,
                    deltax=0.4, lambda = 1e-16){
  Xk = Xstar = matrix(Xinit, ncol=2)
  Cstar = f(Xinit)
  Cvec = rep(0,K)
  Tk = T0
  dk = dmin(Xk)
  cat("dmin:", dk, "\n")
  for(k in 1:K){
    i = tirage(Xk, lambda = lambda) #i = k%%6 +2
    if(k<round(K/10)){
      epsvec = c(-1.2, -0.6,  0, 0.6, 1.2)
    }else{
      if(k<round(K/5)){
        epsvec = c(-1.2, -0.8,-0.4, 0, 0.4, 0.8, 1.2)
      }else{
        epsvec = c(-1.2,-1,-0.8,-0.6, -0.4, -0.2,0,0.2, 0.4,0.6, 0.8,1, 1.2)
      }
    }
    i1 = sample(1:length(epsvec), 1)
    i2 = sample(1:length(epsvec), 1)
    Xk[i,1] = min(max(Xk[i,1] + epsvec[i1], 3), 7)
    Xk[i,2] = min(max(Xk[i,2] + epsvec[i2], 5.4), 14.8)
    dk = dmin(Xk)

    while(dk<deltax){
      i1 = sample(1:length(epsvec), 1)
      i2 = sample(1:length(epsvec), 1)
      Xk[i,1] = min(max(Xk[i,1] + epsvec[i1], 3), 7)
      Xk[i,2] = min(max(Xk[i,2] + epsvec[i2], 5.4),14.8)
      dk = dmin(Xk)
    }
    cat("dmin:", dk, "\n")
    Ck = f(Xk)
    Delta = Cstar - Ck
    if(Delta<0){
      Xstar = Xk
      Cstar = Ck
    }else{
      u = runif(1)
      p = exp(-Delta/Tk)
      if(u<p){
        Xstar = Xk
        Cstar = Ck
      }

    }
    Cvec[k] <- Cstar
    Tk = alpha*Tk
    cat("Iteration RS : ", k, "\n")
  }
  return(list(Xopt=Xstar, Cvec=Cvec))
}


# Plan initial
x1 = c(5, 15) #F1
x8 = c(5, 5.2) #F8
X0 = matrix(0, ncol=2, nrow=8)
X0[1,] = x1
X0[8,] = x8
aa = sample(1:nrow(XX), 6)
X0[2:7,] = XX[aa, ]
maxiter = 1e5
R = RSoptim(f=dmin, Xinit=X0, K=maxiter,
            alpha=0.99, T0 = 1e2, deltax=0)
X0opt = R$Xopt
X0df <- data.frame(
  x =X0opt[,1],
  y =X0opt[,2]
)
ggplot(X0df, aes(x = x, y = y), asp = 1) +
  geom_point(size = 2.5, color = "blue") +
  labs(title = " ", x = "x", y = "y") +
  theme_bw() +
  geom_text(aes(label = ""),
            vjust = -0.5, color = "black") +
  ylim(0, 20) + xlim(0,10)



saveRDS(X0opt, file = "C:/Users/barrya/Desktop/CalibGEOXIM/Datas/X01.rds")

# Plot positions
X0df <- data.frame(
  x = X0opt[,1],
  y = X0opt[,2]
)

png(filename = "X01.png", width = 1920, height = 1920, res = 300)
ggplot(X0df, aes(x = x, y = y), asp = 1) +
  geom_point(size = 2.5, color = "blue") +
  labs(title = " ", x = "x", y = "y") +
  theme_bw() +
  geom_text(aes(label = ""),
            vjust = -0.5, color = "black") +
  ylim(0, 20) + xlim(0,10)
dev.off()

# Evolution optimisation
df <- data.frame(Iteration = 1:length(R$Cvec), Copt = R$Cvec)

png(filename = "EvolRS_dmin.png", width = 1600, height = 1900, res = 300)
ggplot(df, aes(x = Iteration, y = Copt)) +
  geom_line(color = "red", size = 1) +
  labs(title = " ", x = "Itération", y = "dmin") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()


# Optimisation par RS des positions F2,...,F7
DomT = matrix(c(0.2, 0.4, 0.2, -12.31,
                1.6, 4.4, 0.45, -10.39), ncol=2)

x1 = c(5, 15) #F1
x8 = c(5, 5.2) #F8
K=10000
R = RSoptim(f=function(X){CVmodel(X, model, DomT)} ,
            Xinit=X0opt, K=K, alpha=0.99, T0 = 0.01, deltax=1.2)

Xopt = R$Xopt
Cvec = R$Cvec
df <- data.frame(Iteration = 1:K, Copt = R$Cvec[1:K])

# Evolution optimRS
png(filename = "EvolutionC_opt.png", width = 1600, height = 1900, res = 300)
ggplot(df, aes(x = Iteration, y = Copt)) +
  geom_line(color = "blue", size = 1) +
  labs(title = " ", x = "Itération", y = "Copt") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()


# Ordoné le plan
ee = sort.int(Xopt[,2],decreasing = TRUE, index.return = TRUE)$ix
Xopt = Xopt[ee,]
Xopt <- data.frame(
  x =Xopt[,1],
  y =Xopt[,2]
)
xx = seq(0,10, len=50)
yy = seq(0,20, len=100)
Vx <- readRDS("C:/Users/barrya/Desktop/CalibGEOXIM/Datas/Vx.rds")
setwd("C:/Users/barrya/Desktop/CalibGEOXIM/Figures")
png(filename = "Xoptimale.png", width = 1800, height = 1800, res = 175)
image(xx, yy, Vx, col = viridis(100), xlab = "x", ylab = "y",main =" ", asp = 1)
filled.contour(xx, yy, Vx, color.palette = viridis::viridis,
               plot.title = title(main = " ",
                                  xlab = "x", ylab = "y"),
               plot.axes = {
                 axis(1)
                 axis(2)
                 points(Xopt[,1], Xopt[,2], col = "red", pch = 10, lwd = 3)
                 text(Xopt[,1], Xopt[,2], font = 2, cex = 1.75,
                      labels = paste0("F", 1:nrow(Xopt)),
                      pos = 3, col = "white",lwd = 2) })
dev.off()
setwd("C:/Users/barrya/Desktop/CalibGEOXIM/Datas")
saveRDS(Xopt, file = "X_opt.rds")
saveRDS(R, file = "ResultatXopt.rds")

# Observations physiques
iXopt =indx(Xopt, XX)
Yobs <- readRDS("C:/Users/barrya/Desktop/CalibGEOXIM/Datas/YYobs.rds")
Yobs = 1e-7*Yobs
Yobs = Yobs[c(5001:10000), ]
sigeps = 0.1*sqrt(mean(var(Yobs)))
Eps = matrix(rnorm(8*461, mean=0, sd=sigeps), ncol=461)
Yopt = Yobs[iXopt,] + Eps
saveRDS(Yopt, file = "Yopt.rds")
