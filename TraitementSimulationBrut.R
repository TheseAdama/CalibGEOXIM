# Chemin
ds <- "C:/Users/barrya/Desktop/Domain_MolarFraction_NACLw"

# Liste des noms
nmf <- sprintf("Domain_MolarFraction_NACLw%06d", 0:460)

FFsim <- matrix(NA, nrow = 25000, ncol = 461)
for (i in seq_along(nmf)) {
  pathf <- file.path(ds, nmf[i])
  dd <- read.table(pathf, skip = 4, nrows = 25000, header = FALSE)
  FFsim[, i] <- dd[, 1]
  cat("Iteration : ", i, "\n")
}
print(dim(FFsim))
aa = sample(1:25000,100)
tt = seq(0,460, by=1)
plot(1:461, 1e5*FFsim[1,], type='l', lwd=2, ylim=c(11.5,16))
for(j in aa){
  lines(tt, 1e5*FFsim[j,], col=j, lwd=2)
}
head(FFsim)

setwd("C:/Users/barrya/Desktop/CalibGEOXIM/Datas")
saveRDS(FFsim, file = "Fsimthemap.rds")
