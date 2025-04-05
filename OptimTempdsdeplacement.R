# Distribution a posteriori
Pdfl <- function(theta, Xobs, Yobs, sigeps, model, icand){

  # List des indices
  indt = list()
  indt[[1]] = indt[[8]] = c(1:461)
  indt[[2]] = indt[[3]]= indt[[4]] = c(1:icand)
  indt[[5]] = indt[[6]]= indt[[7]] = c((icand+1):461)

  # Prediction Dtheta = (theta,X)
  Dp = cbind(matrix(rep(theta, nrow(Xobs)), ncol=4, byrow = TRUE), Xobs)
  Yp = svdgppredict(model, Dpred=Dp)

  Mm = 1e-5*Yp$Mn
  Mmv = c()
  for(i in 1:length(indt)){
    Mmv <- c(Mmv, as.vector(Mm[i, indt[[i]]]))
  }
  V = matrix(Mmv-Yobs, nrow=1)
  SSm = (1/sigeps)*V%*%t(V)
  P = - 0.5*SSm
  P = exp(as.numeric(P))
  return(P)
}

# Definition du score
Score <- function(icur, icand, Xobs, Yobs, model,
                  Dtheta, Nt=461, K=100, L=100){
  THETA = Unif(L, Dtheta)
  Sig = rgamma(L,shape=4, scale=2.5e-11)
  YY1 = lapply(1:L, function(l){
    Dp = cbind(matrix(THETA[l, ], nrow=1), matrix(Xobs[1,],nrow=1))
    Yp = svdgppredict(model, Dpred=Dp)
    Mm = 1e-5*Yp$Mn
    Eps = matrix(rnorm(ncol(Mm)*nrow(Mm), mean=0, sd=Sig[l]), ncol=ncol(Mm), nrow=nrow(Mm))
    Ys = Mm + Eps
    Ys = as.vector(Ys[,(icur+1):Nt])
    Yobs1 = as.vector(Yobs[1,1:icur])
    R = append(Yobs1,Ys)
    return(R)
  })

  YY8 = lapply(1:L, function(l){
    Dp = cbind(matrix(THETA[l, ], nrow=1), matrix(Xobs[8,],nrow=1))
    Yp = svdgppredict(model, Dpred=Dp)
    Mm = 1e-5*Yp$Mn
    Eps = matrix(rnorm(ncol(Mm)*nrow(Mm), mean=0, sd=Sig[l]), ncol=ncol(Mm), nrow=nrow(Mm))
    Ys = Mm + Eps
    Ys = as.vector(Ys[,(icur+1):Nt])
    Yobs8 = as.vector(Yobs[8, 1:icur])
    R = append(Yobs8,Ys)
    return(R)
  })

  YYG1 = lapply(1:L, function(l){
    if(icur!=icand){
      Dp = cbind(matrix(rep(THETA[l, ], 3), nrow=3, byrow = TRUE), Xobs[2:4,])
      Yp = svdgppredict(model, Dpred=Dp)
      Mm = 1e-5*Yp$Mn
      Eps = matrix(rnorm(ncol(Mm)*nrow(Mm), mean=0, sd=Sig[l]), ncol=ncol(Mm), nrow=nrow(Mm))
      Ys = Mm + Eps
      Ys = Ys[,(icur+1):icand]
      R = matrix(cbind(Yobs[2:4, 1:icur], Ys), nrow=3)
    }else{
      R = matrix(Yobs[2:4, 1:icur], nrow=3)
    }
    R = as.vector(t(R))
    return(R)
  })

  YYG2 = lapply(1:L, function(l){
    Dp = cbind(matrix(rep(THETA[l, ], 3), nrow=3, byrow = TRUE), Xobs[5:7,])
    Yp = svdgppredict(model, Dpred=Dp)
    Mm = 1e-5*Yp$Mn
    Eps = matrix(rnorm(ncol(Mm)*nrow(Mm), mean=0, sd=Sig[l]), ncol=ncol(Mm), nrow=nrow(Mm))
    Ys = Mm + Eps
    Ys = Ys[,(icand+1):Nt]
    Ys = as.vector(t(Ys))
    return(Ys)
  })
  W = E = matrix(0, nrow=L, ncol=K)
  for(l in 1:L){
    TT =  Unif(K, D=Dtheta)
    Ysim = c(YY1[[l]], YYG1[[l]], YYG2[[l]], YY8[[l]])
    W[l, ] = apply(TT, 1, function(theta){
      R = Pdfl(theta, Xobs, Yobs=Ysim, sigeps=Sig[l],model, icand)
      return(R)
    })

    E[l, ] = apply(TT, 1, function(theta){
      V = as.vector(THETA[l, ]-theta)
      R = t(V)%*%V
      return(as.numeric(R)) })

  }

  o = mean(crossprod(W, E))
  return(o)
}


# Optimisation du temps de deplacement
load("C:/Users/barrya/Desktop/Application/model2.RData")
DomT =  matrix(c(0.2, 0.4, 0.2, -12.31,1.6,4.4, 0.45,-10.39), ncol=2)
Xobs <- readRDS("C:/Users/barrya/Desktop/CalibGEOXIM/Xopt.rds")
Yobs <- readRDS("C:/Users/barrya/Desktop/CalibGEOXIM/Yobs.rds")
icurv = seq(23,461, by=23)
L = length(icurv)
count = 1
for(i in 2:L){
  icandt = icurv[i:L]
  Ni = length(icandt)
  ss = rep(0, Ni)
  for(j in 1:Ni){
    ss[j] <- Score(icur=icurv[i], icand=icandt[j], Xobs=Xopt, Yobs=Yopt,
                   model=model2, Dtheta=DomT, Nt=461, K=2, L=2)
    count = count+1
  }
  istar = which.min(ss)
  if(istar==i){
    cat("Temps de changement tg=", icur)
    break
  }
  cat("Iteration : ", i, "\n")

}

cat("Nbre evaluation : ", count, "\n")
saveRDS(icur, file = "topt.rds")
