library(viridis)
library(SVDGP)

load("C:/Users/barrya/Desktop/CalibGEOXIM/Datas/modelm2.RData")
XX <- readRDS("C:/Users/barrya/Desktop/CalibGEOXIM/Datas/XX.rds")
Yobs <- readRDS("C:/Users/barrya/Desktop/CalibGEOXIM/Datas/YYobs.rds")
Yobs = 1e-7*Yobs
Yobs = Yobs[c(5001:10000), ]
Yobs = 1e5*Yobs

#fbar x
Ybar = rep(0,5000)
for(i in 1:5000){
  Ybar[i] = mean(Yobs[i,], na.rm=TRUE)
}

x <- seq(0, 10, len=50)
y <- seq(0, 20, len=100)
z <- matrix(Ybar, ncol=100, nrow=50, byrow=FALSE)
png(filename = "ymean.png", width = 1800, height = 1800, res = 300)
image(x, y, z, col = viridis(100), xlab = "x", ylab = "y",
      main = expression("Moyenne " * bar(f)[x](x) ~ "(10"^-5 ~ " mol/L)"), xaxt = "n", yaxt = "n")
axis(1, at = seq(0, 10, by = 2))
axis(2, at = seq(0, 20, by = 2))
filled.contour(x, y, z, color.palette = viridis,
               plot.title = title(main = expression("Moyenne " * bar(f)[x](x)~ "(10"^-5 ~ " mol/L)"),
                                  xlab = "x", ylab = "y"),
               plot.axes = { axis(1); axis(2) })
dev.off()

# Mubarx
thetamap = c(0.8556, 4.3784, 0.2019,-12.2922)
Fbar = rep(0,5000)
for(i in 1:5000){
  Db = cbind(matrix(thetamap, nrow = 1), matrix(XX[i,], nrow=1))
  YY = svdgppredict(model=model2, Db)
  Fbar[i] = mean(YY$Mn, na.rm=TRUE)
  cat("Iteration",i,"\n")
}

zf <- matrix(Fbar, ncol=100, nrow=50, byrow=FALSE)

png(filename = "mucal.png", width = 1800, height = 1800, res = 300)
image(x, y, zf, col = viridis(100), xlab = "x", ylab = "x",
      main = expression("Moyenne " * bar(mu)[x](x) ~ "(10"^-5 ~ " mol/L)")
      , xaxt = "n", yaxt = "n")
axis(1, at = seq(0, 10, by = 2))
axis(2, at = seq(0, 20, by = 2))
filled.contour(x, y, zf, color.palette = viridis,
               plot.title = title(main =  expression("Moyenne " * bar(mu)[x](x) ~ "(10"^-5 ~ " mol/L)"),               xlab = "x", ylab = "y"),
               plot.axes = { axis(1); axis(2) })
dev.off()

# fbart et mubart
Ybart = rep(0,461)
for(i in 1:461){
  Ybart[i] = mean(Yobs[,i], na.rm=TRUE)
}

FF = matrix(0, ncol=461, nrow=5000)
for(i in 1:5000){
  Db = cbind(matrix(theta0, nrow = 1), matrix(XX[i,], nrow=1))
  YY = svdgppredict(model=model2, Db)
  FF[i, ] = YY$Mn
  cat("Iteration",i,"\n")
}

Ft = rep(0,461)
for(i in 1:461){
  Ft[i] = mean(FF[,i], na.rm=TRUE)
}

data <- data.frame(
  tt = tt,
  Ybart = Ybart,
  Ft = Ft
)
p <- ggplot(data, aes(x = tt)) +
  geom_line(aes(y = Ybart, color = "Moyenne temporelle"), size = 1) +
  geom_line(aes(y = Ft, color = "Prédictions"), size = 1) +
  scale_color_manual(
    name = " ",
    values = c("Moyenne temporelle" = "blue", "Prédictions" = "red"),
    labels = c(
      expression("Moyenne " * bar(f)[t](t)),
      expression("Moyenne " * bar(mu)[t](t))
    )
  )  +
  labs(
    x = expression(Temps ~ (h)),
    y = expression("Valeurs " ~ "(10"^-5 ~ "mol/L)"),
    title = " "
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18),
    legend.position = "top",
    legend.text.align = 0,
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14)
  )

print(p)
ggsave(filename = "ytmut.png", plot = p, width = 10,
       height = 6, dpi = 300)

# Erreur de calibration
theta0 = c(0.869, 4.365, 0.275, -12.25)
thetamap = c(0.8438 , 4.3921 , 0.2538 , -12.3024)
ee = norm(theta0-thetamap, "2")
er = norm((theta0-thetamap), "2")/norm(theta0, "2")
ee
er

sig02 = 6.4566*1e-14
sigmap2 = 3.2981*1e-14
ersig = norm((sig02-sigmap2), "2")/norm(sig02, "2")
ersig
