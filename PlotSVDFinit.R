# Chargement package
library(ggplot2)
library(patchwork)

# Chargement
Dinit <- readRDS("C:/Users/barrya/Desktop/CalibGEOXIM/Dinit.rds")
Finit <- readRDS("C:/Users/barrya/Desktop/CalibGEOXIM/Finit.rds")

dim(Dinit)
dim(Finit)

####--SVD---####
# Valeurs propres spatiaux
tt=seq(0,460, by=1)
A = svd(1e5*Finit)
K = 15
data <- data.frame(
  k = 1:K,
  vp = diag(diag(A$d))[1:K]
)
png(filename = "vp.png", width = 1920, height = 1080, res = 300)
ggplot(data, aes(x = k, y = log(vp))) +
  geom_line(color = "red", size = 1) +
  geom_point(color = "blue", fill = "white", shape = 21, size = 3) +
  labs(title = " ", x = "k", y = expression(log(lambda[k]))) +
  theme_bw()
dev.off()

# Modes temporelles
data <- data.frame(
  Temps = tt,
  V1 = A$v[, 1],
  V2 = A$v[, 2],
  V3 = A$v[, 3],
  V4 = A$v[, 4]
)

# Graphiques 1,2,3,4
plot1 <- ggplot(data, aes(x = Temps, y = V1)) +
  geom_line(color = 'blue', size = 1.2) +
  labs(x = "Temps (h)", y = "V1", title = "Premier mode temporel") +
  theme_bw()

plot2 <- ggplot(data, aes(x = Temps, y = V2)) +
  geom_line(color = 'red', size = 1.2) +
  labs(x = "Temps (h)", y = "V2", title = "Deuxième mode temporel") +
  theme_bw()

plot3 <- ggplot(data, aes(x = Temps, y = V3)) +
  geom_line(color = 'green', size = 1.2) +
  labs(x = "Temps (h)", y = "V3", title = "Troisième mode temporel") +
  theme_bw()

plot4 <- ggplot(data, aes(x = Temps, y = V4)) +
  geom_line(color = 'cyan', size = 1.2) +
  labs(x = "Temps (h)", y = "V4", title = "Quatrième mode temporel") +
  theme_bw()

cbplot <- plot1 + plot2 + plot3 + plot4

png(filename = "Modetemp.png", width = 1920, height = 1080, res = 300)
print(cbplot)
dev.off()

