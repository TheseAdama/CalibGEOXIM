# Bibliothèques
library(viridis)
library(ggplot2)
library(latex2exp)

# Chargement maillage, simulations, plan lhsmaximin theta
XX <- readRDS("C:/Users/barrya/Desktop/CalibGEOXIM/Datas/XX.rds")
FF <- readRDS("C:/Users/barrya/Desktop/CalibGEOXIM/Datas/Array100sim.rds")
TT <- read.table("C:/Users/barrya/Desktop/CalibGEOXIM/Datas/LHS-svg_X_100_saltinj.dat",
                 header = TRUE, sep = " ", quote = "\"", stringsAsFactors = FALSE)
Dth = matrix(c(0.2, 0.4, 0.2,-12.31,1.6,4.4, 0.45,-10.39), ncol=2)
THETA = TT
for(i in 1:dim(TT)[2]){
  THETA[,i] <- (Dth[i,2]-Dth[i,1])*TT[,i] + Dth[i,1]
}

# Selection des simulations
list = c()
for(i in 1:100){
  list[[i]] <- FF[c(5001:10000), ,i]
}
rm(FF)

Fsim = array(0, dim=c(5000,461,100))
for(i in 1:100){
  Fsim[,,i] <- list[[i]]
}
dim(Fsim)

setwd("C:/Users/barrya/Desktop/CalibGEOXIM/Datas")
saveRDS(Fsim, file = "Fsim.rds")
saveRDS(THETA, file = "THETA.rds")

tt = seq(0, 460, by=1)
Dx = matrix(c(0,0,10,20), ncol=2)
Nt = 461
Nx = 5000
Ntheta = 100

# Moyenne et Variation en fonction de t
Moyt = Vart = rep(0, Nt)
for(t in 1:Nt){
  Moyt[t] = mean(Fsim[,t,])
  Vart[t] = sqrt(var(Fsim[,t,]))
  cat("Iteration : ", t, "\n")
}

setwd("C:/Users/barrya/Desktop/CalibGEOXIM/Figures")
## Moyenne
data <- data.frame(x = tt, y = Moyt)
png(filename = "Moyt.png", width = 1920, height = 1080, res = 300)
ggplot(data, aes(x = x, y = 1e5*y)) +
  geom_line(color = "blue2", size = 1) +
  ggtitle("Moyenne Mt") +
  xlab("Temps (h)") +
  ylab(expression(10^{5} ~ Mt)) + #~ "(10"^-5 ~ "mol/L)"
  theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", color = "black"),
    axis.title = element_text(face = "italic"),
    axis.text = element_text(size = 12, color = "black"),
    panel.grid.major = element_line(color = "grey80", linetype = "dotted")
  )
dev.off()

## Variation
data <- data.frame(x = tt, y = Vart)
png(filename = "Vart.png", width = 1920, height = 1080, res = 300)
ggplot(data, aes(x = x, y = 1e5*y)) +
  geom_line(color = "red", size = 1) +
  ggtitle("Variation Vt") +
  xlab("Temps (h)") +
  ylab(expression(10^{5}~Vt)) +
  theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", color = "black"),
    axis.title = element_text(face = "italic"),
    axis.text = element_text(size = 12, color = "black"),
    panel.grid.major = element_line(color = "grey80", linetype = "dotted")
  )
dev.off()

# Quelques simulations tirées aléatoirement
xind = sample(1:5000, 100)
data_list <- list()
for (i in 1:100) {
  ith = sample(1:100)
  data_list[[i]] <- data.frame(
    time = tt,
    F_value = Fsim[xind[i], , sample(1:100,1)],
    sim_id = i
  )
}
aa = sample(1:1000,100)
df <- do.call(rbind, data_list)
palrgb <- colorRampPalette(c("red", "green", "blue", "yellow", "magenta", "cyan"))


png(filename = "Quelquesim.png", width = 1920, height = 1080, res = 300)
ggplot(df, aes(x = time, y = 1e5*F_value, group = sim_id, color = factor(sim_id))) +
  geom_line(size = 1) +
  scale_color_manual(values =palrgb(1000)[aa]) +
  xlab("Temps (h)") +
  ylab(expression("Fraction molaire de NaCl" ~ (x ~ 10^{5}))) + #10^{5} ~ f[code]
  ylim(11.4, 13.75) +
  theme_bw() +
  guides(color = "none")
dev.off()


# Moyenne et Variation en fonction des positions
Moyx = Varx = rep(0, Nx)
for(i in 1:Nx){
  Moyx[i] = mean(Fsim[i, , ], na.rm = TRUE)
  Varx[i] = sqrt(var(Fsim[i, , ],na.rm = TRUE))
  cat("Iteration : ", i, "\n")
}
Moyx = matrix(Moyx, nrow=50, ncol=100, byrow = FALSE)
Varx = matrix(Varx, nrow=50, ncol=100, byrow = FALSE)
x = seq(0,10, length=50)
y = seq(0,20, length=100)

## Moyenne
png(filename = "Moyx.png", width = 1800, height = 1900, res = 300)
image(x, y, 1e5*Moyx, col = viridis(100),
      xlab = "X (m)", ylab = "Y (m)",
      main =expression(10^{5}~Mx), asp = 1)
filled.contour(x, y, 1e5*Moyx, color.palette = viridis,
               plot.title = title(main = expression(10^{5}~Mx),
                                  xlab = "x", ylab = "y"),
               plot.axes = {axis(1); axis(2)},
               key.title = expression(10^{5}~Mx))
dev.off()

## Variation
png(filename = "Varx.png", width = 1800, height = 1900, res = 300)
image(x, y, 1e5*Varx, col = viridis(100),
      xlab = "x (m)", ylab = "y (m)",
      main =expression(10^{5}~Vx), asp = 1)
filled.contour(x, y, 1e5*Varx, color.palette = viridis,
               plot.title = title(main = expression(10^{5}~Vx),
                                  xlab = "x", ylab = "y"),
               plot.axes = {axis(1); axis(2)},
               key.title = expression(10^{5}~Vx))
dev.off()
