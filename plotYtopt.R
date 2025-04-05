library(viridis)
library(ggplot2)
library(dplyr)

Yopt <- readRDS("C:/Users/barrya/Desktop/CalibGEOXIM/Datas/Yopt.rds")
topt <- readRDS("C:/Users/barrya/Desktop/CalibGEOXIM/Datas/topt.rds")
tt <- seq(0, 460, len = 461)
tt1 <- seq(0, topt, by = 1)
tt2 <- seq(topt + 1, 460, by = 1)

Yobs1 <- Yopt[1,]
Yobs2 <- Yopt[2, 1:(topt + 1)]
Yobs3 <- Yopt[3, 1:(topt + 1)]
Yobs4 <- Yopt[4, 1:(topt + 1)]
Yobs5 <- Yopt[5, (topt + 2):461]
Yobs6 <- Yopt[6, (topt + 2):461]
Yobs7 <- Yopt[7, (topt + 2):461]
Yobs8 <- Yopt[8,]

data1 <- data.frame(
  time = c(tt, tt1, tt1, tt1),
  value = c(Yobs1, Yobs2, Yobs3, Yobs4),
  group = factor(rep(1:4, times = c(length(tt), length(tt1), length(tt1), length(tt1))))
)
setwd("C:/Users/barrya/Desktop/CalibGEOXIM/Figures")
png(filename = "Ytopt1.png", width = 1920, height = 1080, res = 300)
ggplot(data1, aes(x = time, y = 1e5*value, color = group)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("blue", "red", "green", "purple")) +
  labs(x = "Temps (h)",
       y = expression("Fraction molaire de NaCl" ~ (x ~ 10^{5})), color = "Puits") +
  xlim(0, 460) +
  theme_bw()
dev.off()

data2 <- data.frame(
  time = c(tt2, tt2, tt2, tt),
  value = c(Yobs5, Yobs6, Yobs7, Yobs8),
  group = factor(rep(5:8, times = c(length(tt2), length(tt2), length(tt2), length(tt))))
)


png(filename = "Ytopt2.png", width = 1920, height = 1080, res = 300)
ggplot(data2, aes(x = time, y = 1e5*value, color = group)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = c("orange", "brown", "pink", "cyan")) +
  labs(x = "Temps (h)", y = expression("Fraction molaire de NaCl" ~ (x ~ 10^{5})), color = "Puits") +
  xlim(0, 460) +
  theme_bw()
dev.off()

