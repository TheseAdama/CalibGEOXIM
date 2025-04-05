library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)

# Chargement modele
load("C:/Users/barrya/Desktop/CalibGEOXIM/Datas/modelm0.RData")
setwd("C:/Users/barrya/Desktop/CalibGEOXIM/Figures")

# Prediction
Dpred = Dxtr[sample(1:nrow(Dxtr), 1000),]
FDpred = 1e5*SimCorresp(Dpred, XX, THETA, Fsim, Nt=461)
YY = svdgppredict(model=model2, Dpred)
Mn = YY$Mn
Qg2 = 1 - (sum((FDpred - Mn)^ 2)/sum((FDpred - mean(FDpred))^2))
print(100*Qg2)

# Evaluation des sous modeles
U = (FDpred-colMeans(FDpred))%*%model2$V
Q2 = rep(0, K)
for(k in 1:K){
  Q2[k] = 1 - (sum((U[,k] - YY$Ypred[,k])^ 2)/sum((U[,k] - mean(U[,k]))^2))
}
print(100*Q2)

# modele GP1
data <- data.frame(X = U[,1], Y = YY$Ypred[,1])
png(filename = "GadeqU1.png", width = 1200, height = 1200, res = 300)
ggplot(data, aes(x = X, y = Y)) +
  geom_point(color = "blue", size = 1) +
  geom_abline(slope = 1, intercept = 0, color = "red", size = 1.25) + # Trait plein sans linetype
  labs(x = "Vraies valeurs", y = "Valeurs prédites",
       title = expression("Graphique d'adéquation GP1 " ~ Q[2] == 95.16 *"%")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

# Modele GP2
data <- data.frame(X = U[,2], Y = YY$Ypred[,2])
png(filename = "GadeqU2.png", width = 1200, height = 1200, res = 300)
ggplot(data, aes(x = X, y = Y)) +
  geom_point(color = "blue", size = 1.5) +
  geom_abline(slope = 1, intercept = 0, color = "red", size = 1.25) + # Trait plein sans linetype
  labs(x = "Vraies valeurs", y = "Valeurs prédites",
       title = expression("Graphique d'adéquation GP2 " ~ Q[2] == 88.48 *"%")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

# Modele GP3
data <- data.frame(X = U[,3], Y = YY$Ypred[,3])
png(filename = "GadeqU3.png", width = 1200, height = 1200, res = 300)
ggplot(data, aes(x = X, y = Y)) +
  geom_point(color = "blue", size = 1.5) +
  geom_abline(slope = 1, intercept = 0, color = "red", size = 1.25) + # Trait plein sans linetype
  labs(x = "Vraies valeurs", y = "Valeurs prédites",
       title = expression("Graphique d'adéquation GP3 " ~ Q[2] == 72.85 *"%")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

# Modele GP4
data <- data.frame(X = U[,4], Y = YY$Ypred[,4])

png(filename = "GadeqU4.png", width = 1200, height = 1200, res = 300)
ggplot(data, aes(x = X, y = Y)) +
  geom_point(color = "blue", size = 1.5) +
  geom_abline(slope = 1, intercept = 0, color = "red", size = 1.25) + # Trait plein sans linetype
  labs(x = "Vraies valeurs", y = "Valeurs prédites",
       title = expression("Graphique d'adéquation GP4 " ~ Q[2] == 66.85 *"%")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

# Graphique des différences VRAIE- PREDICTION
# Échantillonnage
Df <- (FDpred - Mn)/FDpred
set.seed(123)
aa <- sample(1:nrow(Df), 100)
Df <- Df[aa, ]

# Conversion des données en un format long pour ggplot
df_long <- as.data.frame(Df) %>%
  mutate(ID = row_number()) %>%
  pivot_longer(-ID, names_to = "Time", values_to = "Difference") %>%
  mutate(Time = as.numeric(gsub("V", "", Time)),
         Time = tt[Time])

# Graphique
png(filename = "Diffadeq.png", width = 1920, height = 1200, res = 300)
ggplot(df_long, aes(x = Time, y = Difference, group = ID, color = as.factor(ID))) +
  geom_line(size = 1) +
  scale_color_viridis_d(option = "plasma", guide = FALSE) +
  labs(x = "Temps (h)", y = expression("Différence relative"),
       title = " ") +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10)) +
  ylim(-0.0665, 0.0225)
dev.off()

# Graphique d'adéquation
FF <- as.vector(colMeans(FDpred))
YY <- as.vector(colMeans(Mn))
data <- data.frame(FF = FF, YY = YY)

png(filename = "Gadeq.png", width = 1920, height = 1100, res = 300)
ggplot(data, aes(x = FF, y = YY)) +
  geom_point(color = "blue", size = 1.5) +
  geom_abline(slope = 1, intercept = 0, color = "red", size = 1.5) +
  labs(x = "Simulations", y = "Prédictions",
       title = "Graphique d'adéquation (Q2 = 82.43%)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()


# Q2(x,y)
# Charger le modele
load("C:/Users/barrya/Desktop/CalibGEOXIM/Datas/model0.RData")
Nx = 5000
Nt = 461
Qx2 =rep(0, Nx)
for(i in 1:Nx){
  x = XX[i, ]
  Dpred = cbind(THETA, matrix(rep(x, 100), ncol=2, byrow = TRUE))
  FDpred = 1e5*SimCorresp(Dpred, XX, THETA, Fsim, Nt=461)
  YY = svdgppredict(model, Dpred)
  Mn = YY$Mn
  Qx2[i] = 1 - (sum((FDpred - Mn)^ 2)/sum((FDpred - colMeans(FDpred))^2))
  cat("Iteration : ", i, "\n")
}
Qx2 = matrix(Qx2, nrow = 50, ncol = 100, byrow = FALSE)
x = seq(0, 10, len=50)
y = seq(0, 20, len=100)

png(filename = "Q2m1.png", width = 1920, height = 1080, res = 300)
image(x, y, Qx2, col = viridis(100), xlab = "x", ylab = "y", main ="Q2", asp = 1)
filled.contour(x, y, Qx2, color.palette = viridis,
               plot.title = title(main ="Q2", xlab = "x", ylab = "y"),
               plot.axes = { axis(1); axis(2) })

setwd("C:/Users/barrya/Desktop/CalibGEOXIM/Datas")
saveRDS(Qx2, "Q2m1.rds")

# Charger le modele 2
load("C:/Users/barrya/Desktop/CalibGEOXIM/Datas/modelm2.RData")
# Calcul Q2 spatial
Nx = 5000
Nt = 461
Qx2 =rep(0, Nx)
for(i in 1:Nx){
  x = XX[i, ]
  Dpred = cbind(THETA, matrix(rep(x, 100), ncol=2, byrow = TRUE))
  FDpred = 1e4*SimCorresp(Dpred, XX, THETA, Fsim, Nt=461)
  YY = svdgppredict(model2, Dpred)
  Mn = YY$Mn
  Qx2[i] = 1 - (sum((FDpred - Mn)^ 2)/sum((FDpred - colMeans(FDpred))^2))
  cat("Iteration : ", i, "\n")
}
Qx2 = matrix(Qx2, nrow = 50, ncol = 100, byrow = FALSE)
x = seq(0, 10, len=50)
y = seq(0, 20, len=100)

png(filename = "Q2m2.png", width = 1800, height = 1900, res = 300)
image(x, y, Qx2, col = viridis(100), xlab = "x", ylab = "y", main ="Q2", asp = 1)
filled.contour(x, y, Qx2, color.palette = viridis,
               plot.title = title(main ="Q2", xlab = "x", ylab = "y"),
               plot.axes = { axis(1); axis(2) })
dev.off()

setwd("C:/Users/barrya/Desktop/CalibGEOXIM/Datas")
saveRDS(Qx2, "Q2m2.rds")
