library(coda)
Rn <- readRDS("C:/Users/barrya/Desktop/CalibGEOXIM/Datas/sampleMCMC3.rds")
namepar <- c("theta1", "theta2","theta3","theta4","sigma2" )
as <- 4.56e-16
bs <- 1.82e-13
D <- matrix(c(0.2, 0.4, 0.2, -12.31, 4.56e-16,
              1.6, 4.4, 0.45, -10.39, 1.82e-13), ncol = 2)

setwd("C:/Users/barrya/Desktop/CalibGEOXIM/Figures")
png(filename = "ACF1.png", width = 1920, height = 1920, res = 300)
par(mfrow = c(2, 3))
for(i in 1:5){
  acf(Rn[, i], main = bquote("ACF pour" ~ .(namepar[i])))
}
dev.off()

png(filename = "ACF2.png", width = 1920, height = 1920, res = 300)
par(mfrow = c(2, 3))
for(i in 1:5){
  acf(Rn[seq(1, nrow(Rn), by=30), i], main = paste("ACF pour", namepar[i]))
}
dev.off()

R = Rn[seq(1, nrow(Rn), by=30),]

png(filename = "Thetachain1.png", width = 1920, height = 1920, res = 300)
par(mfrow = c(2, 2))
for(i in 1:4) {
  plot(R[, i], type = "l", main = "Chaîne MCMC ",
       ylim = c(min(R[,i])*0.99, max(R[,i])*1.01),
       ylab = bquote(theta[.(i)]), xlab = "Itérations")
}
dev.off()

png(filename = "varepschain.png", width = 1200, height = 1200, res = 300)
par(mfrow = c(1, 1))
plot(R[,5], type = "l", main = "Chaîne MCMC",
     ylim = c(min(R[,5])*0.98, max(R[,5])*1.02),
     ylab = expression(sigma[epsilon]^2), xlab = "Itérations")
dev.off()

rs <- function(R,i){
  dens = density(R[, i], bw = 2*density(R[, i])$bw, from = D[i, 1], to = D[i, 2])
  dd <- data.frame(mean=mean(R[,i]),
                   var = var(R[,i]),
                   med = median(R[,i]),
                   maxp = dens$x[which.max(dens$y)],
                   qmin = quantile(R[,i], 0.025),
                   qmax = quantile(R[,i], 0.975)
  )
  return(dd)
}

d1 = rs(Rn,1)
d2 = rs(Rn,2)
d3 = rs(Rn,3)
d4 = rs(Rn,4)
d5 = rs(Rn,5)
d1
d2
d3
d4
d5

