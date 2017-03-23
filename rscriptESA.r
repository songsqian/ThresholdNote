## Outline:
## 1. Methods in Qian et al (2003) are computational methods for solving a specific model
## 2. Any model is an assumption on the likely underlying data generating process
## 3. Whether a model fits a specific data must be verified, often based on the goodness-of-fit of the fitted model to data
##    and the specific probabilistic assumptions
## 4. Statistical models can be evaluated by examing the residuals
## 5. When a model is selected, model parameters are estimated such that the resulting model's fit to the data is optimal.
##    As a result, it is often difficult to judge the appropriateness of a model 
## 6. Alternative models must be evaluated: add a linear model

## Computational details:
### bugs input:

ini <- function(mat, val) {                         # Initializing matrices
  ina <- is.na(mat); mat[ina] <- val; mat[!ina] <- NA; mat;}

## the general model
writemodel <- 1
if (writemodel){
  cat("
## the general model
model{
  for (i in 1:n){
    y[i] ~ dnorm(y.hat[i], prec)
    y.hat[i] <- (beta[1]+delta[1]*step(x[i]-phi)) +
                (beta[2]+delta[2]*step(x[i]-phi))*(x[i]-phi)
  }
  for (i in 1:2){
    delta[i] ~ dnorm(0, 0.01)
    beta[i] ~ dnorm(0, 0.001)
  }
  phi ~ dunif(l,u)
  prec <- pow(sigma, -2)
  sigma ~ dunif(0, 10)
}
  ", file="threshG.txt");
}

bugs.in <- function(infile, y.col, x.col,
                    y.trans="none", x.trans="none", Subset=NULL){
  if(!is.null(Subset)) infile<-infile[Subset,]
  if (y.trans=="logit") y <- logit(infile[,y.col])
  else if (y.trans=="log") y <- log(infile[,y.col])
  else y <- infile[,y.col]
  n <- length(y)
  if (x.trans=="logit") x <- logit(infile[,x.col])
  else if (x.trans=="log") x <- log(infile[,x.col])
  else x <- infile[,x.col]
  bugs.dat <- list(n=n, y=y, x=x, l=min(x), u=max(x))
  init1 <- list(delta=rnorm(2),
                beta=rnorm(2), phi=runif(1, min(x),max(x)),
                sigma=runif(1), y=ini(y, rnorm(1, 6, 1.4)))
  init2 <- list(delta=rnorm(2),
                beta=rnorm(2), phi=runif(1, min(x),max(x)),
                sigma=runif(1), y=ini(y, rnorm(1, 6, 1.4)))
  init3 <- list(delta=rnorm(2),
                beta=rnorm(2), phi=runif(1, min(x),max(x)),
                sigma=runif(1), y=ini(y, rnorm(1, 6, 1.4)))
  inits <- list(init1, init2, init3)
  para <- c("beta","delta","phi", "sigma")
  return(list(data=bugs.dat, inits=inits, para=para))
}
#
##


## for the hockey stick model
if (writemodel){
  cat("## hockey stick
model{
  for (i in 1:n){
    y[i] ~ dnorm(y.hat[i], prec)
    y.hat[i] <- beta[1] +
             (beta[2]+delta*step(x[i]-phi))*(x[i]-phi)
  }
  for (i in 1:2){
    beta[i] ~ dnorm(0, 0.0001)
  }
  delta ~ dnorm(0, 0.0001)
  phi ~ dunif(l,u)
  prec <- pow(sigma, -2)
  sigma ~ dunif(0, 10)
}
  ", file="threshH.txt");
}

bugs.in2 <- function(infile, y.col, x.col,
                     y.trans="none", x.trans="none", Subset=NULL){
  if(!is.null(Subset)) infile<-infile[Subset,]
  if (y.trans=="logit") y <- logit(infile[,y.col])
  else if (y.trans=="log") y <- log(infile[,y.col])
  else y <- infile[,y.col]
  n <- length(y)
  if (x.trans=="logit") x <- logit(infile[,x.col])
  else if (x.trans=="log") x <- log(infile[,x.col])
  else x <- infile[,x.col]
  bugs.dat <- list(n=n, y=y, x=x, l=min(x), u=max(x))
  init1 <- list(delta=rnorm(1),
                beta=rnorm(2), phi=runif(1, min(x),max(x)),
                sigma=runif(1), y=ini(y, rnorm(1, 6, 1.4)))
  init2 <- list(delta=rnorm(1),
                beta=rnorm(2), phi=runif(1, min(x),max(x)),
                sigma=runif(1), y=ini(y, rnorm(1, 6, 1.4)))
  init3 <- list(delta=rnorm(1),
                beta=rnorm(2), phi=runif(1, min(x),max(x)),
                sigma=runif(1), y=ini(y, rnorm(1, 6, 1.4)))
  inits <- list(init1, init2, init3)
  para <- c("beta","delta","phi", "sigma")
  return(list(data=bugs.dat, inits=inits, para=para))
}
#
##

### for step function
if (writemodel){
  cat("## the step function
model{
  for (i in 1:n){
    y[i] ~ dnorm(y.hat[i], prec[J[i]])
      y.hat[i] <- beta +
             delta*step(x[i]-phi)
    J[i] <- 1+step(x[i]-phi)
  }
    delta ~ dnorm(0, 0.0001)
    beta ~ dnorm(0, 0.0001)
  phi ~ dunif(l,u)
  for (j in 1:2){
    prec[j] <- pow(sigma[j], -2)
    sigma[j] ~ dunif(0, 50)
  }
}
  ", file="threshS.txt");
}

bugs.in3 <- function(infile, y.col, x.col,
                     y.trans="none", x.trans="none",Subset=NULL){
  if(!is.null(Subset)) infile<-infile[Subset,]
  if (y.trans=="logit") y <- logit(infile[,y.col])
  else if (y.trans=="log") y <- log(infile[,y.col])
  else y <- infile[,y.col]
  n <- length(y)
  if (x.trans=="logit") x <- logit(infile[,x.col])
  else if (x.trans=="log") x <- log(infile[,x.col])
  else x <- infile[,x.col]
  bugs.dat <- list(n=n, y=y, x=x, l=min(x, na.rm=T), u=max(x, na.rm=T))
  init1 <- list(delta=rnorm(1),
                beta=rnorm(1), phi=runif(1, min(x),max(x)),
                sigma=runif(2), y=ini(y, rnorm(1, 6, 1.4)))
  init2 <- list(delta=rnorm(1),
                beta=rnorm(1), phi=runif(1, min(x),max(x)),
                sigma=runif(2), y=ini(y, rnorm(1, 6, 1.4)))
  init3 <- list(delta=rnorm(1),
                beta=rnorm(1), phi=runif(1, min(x),max(x)),
                sigma=runif(2), y=ini(y, rnorm(1, 6, 1.4)))
  inits <- list(init1, init2, init3)
  para <- c("beta","delta","phi", "sigma")
  return(list(data=bugs.dat, inits=inits, para=para))
}

# for linear model
writemodel <- 1
if (writemodel){
  cat("## linear model
  model{
    for (i in 1:n){
      y[i] ~ dnorm(mu[i], tau);
      mu[i] <- beta[1] + beta[2]*x[i];
    }
    for (j in 1:2){
      beta[j] ~ dnorm(0, 0.0001);
    }
    tau <- pow(sigma, -2);
    sigma ~ dunif(0, 10);
  }
  ", file="threshL.txt");
}


bugs.in4 <- function(infile, y.col, x.col,
                     y.trans="none", x.trans="none",Subset=NULL){
  if(!is.null(Subset)) infile<-infile[Subset,]
  if (y.trans=="logit") y <- logit(infile[,y.col])
  else if (y.trans=="log") y <- log(infile[,y.col])
  else y <- infile[,y.col]
  n <- length(y)
  if (x.trans=="logit") x <- logit(infile[,x.col])
  else if (x.trans=="log") x <- log(infile[,x.col])
  else x <- infile[,x.col]
  bugs.dat <- list(n=n, y=y, x=x)
  inits<-list()
  for (i in 1:3)
    inits[[i]] <- list(beta=rnorm(2), sigma=runif(1),
                       y=ini(y, rnorm(1, mean(y,na.rm=T), sd(y,na.rm=T))))
  para <- c("beta","sigma")
  return(list(data=bugs.dat, inits=inits, para=para))
}


## a list of model runs -- checking for threshold responses
##   the prior distribution for the threshold is uniform within
##   the predictor data range.  if the posterior distribution is
##   bimodal at both ends, there is no threshold.


model.runs <- function(infile, y.col, x.col,
                       y.trans="none", x.trans="none", Subset=NULL,
                       xlab="", main="", n.iters=25000, n.keep=1000, n.chains=3){
  print("+++++++++++++++++++++++++++++++++++++++++")
  print(paste("y = ", y.col, "; ", "x = ", x.col))
  infile <- infile[!is.na(infile[,x.col]), ]
  input.to.bugs <- bugs.in(infile=infile,y.col=y.col, x.col=x.col,
                           Subset=Subset, x.trans=x.trans, y.trans=y.trans)
#  print(input.to.bugs$data)
  m <- jags.model("threshG.txt",
                  input.to.bugs$data,
                  input.to.bugs$inits,
                  n.chains=n.chains)
  update(m, n.iters)
  x <- coda.samples(m, input.to.bugs$para,
                    n.iter=n.iters, thin=max(c(1, floor(n.iters*n.chains/n.keep))))
  simCoda1 <- NULL
  for (i in 1:n.chains)
    simCoda1 <- rbind(simCoda1, x[[i]])
 # print(summary(simCoda1))
  ycol <- labels(simCoda1)[[2]]=="phi"
  hist(simCoda1[,ycol], xlab=xlab, main=paste(main, "(general)"),
       xlim=range(input.to.bugs$data$x, na.rm=T), col="gray",
       axes=F, prob=T)
  axis(1)
#  lines(density(simCoda1[,3], from=min(input.to.bugs$data$x),
#               to=max(input.to.bugs$data$x), bw = "sj"))

### ii. the hockey stick model

  input.to.bugs <- bugs.in2(infile=infile,y.col=y.col, x.col=x.col,
                            Subset=Subset, x.trans=x.trans, y.trans=y.trans)
  m <- jags.model("threshH.txt",
                  input.to.bugs$data,
                  input.to.bugs$inits,
                  n.chains=n.chains)
  update(m, n.iters)
  x <- coda.samples(m, input.to.bugs$para,
                    n.iter=n.iters, thin=max(c(1, floor(n.iters*n.chains/n.keep))))

  simCoda2 <- NULL
  for (i in 1:n.chains)
    simCoda2 <- rbind(simCoda2, x[[i]])

  ycol <- labels(simCoda2)[[2]]=="phi"
  hist(simCoda2[,ycol], xlab=xlab, main=paste(main, "(hockey)"),
      xlim=range(input.to.bugs$data$x, na.rm=T), col="gray",
      axes=F, prob=T)
  axis(1)

### iii the step function model

  input.to.bugs <- bugs.in3(infile=infile,y.col=y.col, x.col=x.col,
                            Subset=Subset, x.trans=x.trans, y.trans=y.trans)
  m <- jags.model("threshS.txt",
                  input.to.bugs$data,
                  input.to.bugs$inits,
                  n.chains=n.chains)
  update(m, n.iters)
  x <- coda.samples(m, input.to.bugs$para,
                    n.iter=n.iters, thin=max(c(1, floor(n.iters*n.chains/n.keep))))

  simCoda3 <- NULL
  for (i in 1:n.chains)
    simCoda3 <- rbind(simCoda3, x[[i]])

  ycol<-labels(simCoda3)[[2]]=="phi"
  hist(simCoda3[,ycol], xlab=xlab, main=paste(main, "(step)"),
       xlim=range(input.to.bugs$data$x, na.rm=T), col="gray",
       axes=F, prob=T)
  axis(1)

### iii the step function model

  input.to.bugs <- bugs.in4(infile=infile,y.col=y.col, x.col=x.col,
                            Subset=Subset, x.trans=x.trans, y.trans=y.trans)
  m <- jags.model("threshL.txt",
                  input.to.bugs$data,
                  input.to.bugs$inits,
                  n.chains=n.chains)
  update(m, n.iters)
  x <- coda.samples(m, input.to.bugs$para,
                    n.iter=n.iters, thin=max(c(1, floor(n.iters*n.chains/n.keep))))

  simCoda4 <- NULL
  for (i in 1:n.chains)
    simCoda4 <- rbind(simCoda4, x[[i]])

  invisible(list("gen"=simCoda1, "hockey"=simCoda2, "step"=simCoda3,
                 "linear"=simCoda4))
}


### required packages ###
require(rjags)
require(rv)
require(lattice)
 source("~/Dropbox/rutil/util.R")
 base <- "~/SpiderOak/Threshold/EcolMod2"

## Qian et al 2003 used BCD and % tolerance species

setwd(base)
dataDIROhio <- paste(base, "data", sep="/")

plotDIR <- paste(base, "manuscript", sep="/")

dataDIR <- paste(base, "Data", sep="/")
diatoms <- read.csv(paste(dataDIR, "DiatomGeomeanUTP2mo.csv", sep="/"),
                    header=T, na.string="na")
macroinv <- read.csv(paste(dataDIR, "InvertebratesGeomeansUTP6mo.csv", sep="/"),
                     header=T, na.string="na")

macroinv.tolerant <- read.csv(paste(dataDIR, "macroinv.tolerant.csv",sep="/"),
                              sep=",", header=T)
macroinv.tolerant$X.tolerant <- macroinv.tolerant$X.tolerant/100
## Diatoms
## 1. separating substrates
xyplot(Percent ~ GM.UTP|HABITAT, data=diatoms, subset=GM.UTP<300,
	   xlab="TP (ppb)", ylab="% diatoms")

xyplot(logit(Percent) ~ log(GM.UTP)|HABITAT, data=diatoms, subset=GM.UTP<300,
	   xlab="log TP (ppb)", ylab="logit % diatoms")

subS <- levels(diatoms$HABITAT)
temp1 <- list()
par(mfrow=c(1,3), mar=c(3,0.75,3,0.75), mgp=c(1.5,.5,0))
for (i in 1:length(subS))
temp1[[i]] <- model.runs (infile=diatoms[diatoms$HABITAT==subS[i] &
                            (!is.na(diatoms$GM.UTP)),], y.col="Percent", x.col="GM.UTP",
            y.trans="logit", x.trans="log",
            xlab="x", main=subS, n.iters=100000)

dput(temp1, file="diatomJAGSout.Rdata")

xyplot(BCD ~ GM.UTP|DATE, data=macroinv, subset=GM.UTP<300,
	   xlab="TP (ppb)", ylab="BCD")

xyplot(BCD ~ log(GM.UTP)|DATE, data=macroinv, subset=GM.UTP<300,
	   xlab="log TP (ppb)", ylab="BCD")

## BCD
bcdD <- levels(macroinv$DATE)
temp2 <- list()
par(mfrow=c(1,3), mar=c(3,0.75,3,0.75), mgp=c(1.5,.5,0))
for (i in 1:length(bcdD))
temp2[[i]] <- model.runs (infile=macroinv[macroinv$DATE==bcdD[i] &
                            (!is.na(macroinv$GM.UTP)),], y.col="BCD", x.col="GM.UTP",
            y.trans="none", x.trans="log",
            xlab="x", main=bcdD[i], n.iters=100000)

dput(temp2, file="bcdJAGSout.Rdata")

## DEP SENS
xyplot(DEP.SENS ~(GM.UTP)|DATE, data=macroinv, subset=GM.UTP<300,
	   xlab="TP (ppb)", ylab="% Sensitive species")

xyplot(logit(DEP.SENS) ~ log(GM.UTP)|DATE, data=macroinv, subset=GM.UTP<300,
	   xlab="log TP (ppb)", ylab="logit % sensitive")

temp3 <- list()
par(mfrow=c(1,3), mar=c(3,0.75,3,0.75), mgp=c(1.5,.5,0))
for (i in 1:length(bcdD))
temp3[[i]] <- model.runs (infile=macroinv[macroinv$DATE==bcdD[i] &
                            (!is.na(macroinv$GM.UTP)),], y.col="DEP.SENS", x.col="GM.UTP",
            y.trans="logit", x.trans="log",
            xlab="x", main=bcdD[i], n.iters=100000)

dput(temp3, file="depsensJAGSout.Rdata")

temp4 <- list()
par(mfrow=c(1,3), mar=c(3,0.75,3,0.75), mgp=c(1.5,.5,0))
for (i in 1:length(bcdD))
temp4[[i]] <- model.runs (infile=macroinv[macroinv$DATE==bcdD[i] &
                            (!is.na(macroinv$GM.UTP)),], y.col="DEP.SENS", x.col="GM.UTP",
            y.trans="none", x.trans="log",
            xlab="x", main=bcdD[i], n.iters=100000)

dput(temp4, file="depsensJAGSout2.Rdata")

## tolerant species

tolD <- levels(macroinv.tolerant$Date)[c(3,1,2,4)]
temp5 <- list()
par(mfrow=c(1,3), mar=c(3,0.75,3,0.75), mgp=c(1.5,.5,0))
for (i in 1:length(tolD))
temp5[[i]] <- model.runs (infile=macroinv.tolerant[macroinv.tolerant$Date==tolD[i],],
                          y.col="X.tolerant", x.col="TP_6mo_gm",
                          y.trans="none", x.trans="log",
                          xlab="log TP", main=tolD[i], n.iters=100000)

dput(temp5, file="tolerantJAGSout.Rdata")

temp6 <- list()
par(mfrow=c(1,3), mar=c(3,0.75,3,0.75), mgp=c(1.5,.5,0))
for (i in 1:length(tolD))
temp6[[i]] <- model.runs (infile=macroinv.tolerant[macroinv.tolerant$Date==tolD[i],],
                          y.col="X.tolerant", x.col="TP_6mo_gm",
                          y.trans="logit", x.trans="log",
                          xlab="log TP", main=tolD[i], n.iters=100000)

dput(temp6, file="tolerantJAGSout2.Rdata")


dget (paste(base, "diatomJAGSout.Rdata", sep="/"))->simdataDiatom
dget (paste(base, "bcdJAGSout.Rdata", sep="/"))->simdataBCD
dget (paste(base, "depsensJAGSout.Rdata", sep="/"))->simdataSENS
dget (paste(base, "depsensJAGSout2.Rdata", sep="/"))->simdataSENS2
dget (paste(base, "tolerantJAGSout.Rdata", sep="/"))->simdataTOL
dget (paste(base, "tolerantJAGSout2.Rdata", sep="/"))->simdataTOL2



## plot residuals

## residual plots and Bayesian p-values
## general model:
##     y.hat[i] <- (beta[1]+delta[1]*step(x[i]-phi)) +
##                 (beta[2]+delta[2]*step(x[i]-phi))*(x[i]-phi)
##
## hockey stick model:
##     y.hat[i] <- beta[1] + (beta[2]+delta*step(x[i]-phi))*(x[i]-phi)
##

## step function model
## y[i] ~ dnorm(y.hat[i], prec[J[i]]) y.hat[i] <- beta + delta*step(x[i]-phi)

## 1. Diatom -- subS[1]

## create a plotting function
threshold.plots <- function(obsdata=diatoms[diatoms$HABITAT==subS[2] &
                  (!is.na(diatoms$GM.UTP)),],Txt1="Diatoms", Txt2=subS[2],
                y.col="Percent", x.col="GM.UTP", y.trans=logit, x.trans=log,
                simdata=simdataDiatom[[2]], Xlab="log TP", Ylab="logit % diatom"){
  ## extracting model results
  model1 <- rvsims(simdata[[1]]) ## gen
  model2 <- rvsims(simdata[[2]]) ## hockey
  model3 <- rvsims(simdata[[3]]) ## step
  model4 <- rvsims(simdata[[4]]) ## linear
print(1)
  ## extracting data
  y1 <- y.trans(obsdata[,y.col])
  x1 <- x.trans(obsdata[,x.col])
  resids1<-y1-(yhat1<-(model1["beta[1]"]+model1["delta[1]"]*(x1>=model1["phi"]))+
               (model1["beta[2]"]+model1["delta[2]"]*(x1>=model1["phi"]))*
               (x1-model1["phi"]))
print(2)
  resids2<-y1-(yhat2<-model2["beta[1]"] +
               (model2["beta[2]"]+model2["delta"]*(x1>=model2["phi"]))*
               (x1-model2["phi"]))
print(3)
  resids3<-y1-(yhat3<-model3["beta"] + model3["delta"]*(x1>=model3["phi"]))
print(4)
  resids4<-y1-(yhat4<-model4["beta[1]"] + model4["beta[2]"]*x1)
print(5)
  ## plots
  par(mfrow=c(3,4), oma=c(3, 3, 2, 0.25), mar=c(0.25, 0.25, 0.25, 0.25), mgp=c(1.25,0.125,0.),
      las=1, tck=0.01)

  #### residuals
  sumR1 <- summary(resids1)
  sumR2 <- summary(resids2)
  sumR3 <- summary(resids3)
  sumR4 <- summary(resids4)

  ylims <- range(sumR1[,c(4,8)], sumR2[,c(4,8)], sumR3[,c(4,8)], sumR4[,c(4,8)])
  plot(range(x1), ylims, type="n", xlab="", ylab="", axes=F)
  segments(x0=x1, x1=x1, y0=sumR1[,4], y1=sumR1[,8])
  segments(x0=x1, x1=x1, y0=sumR1[,5], y1=sumR1[,7], lwd=3)
  points(x1, sumR1[,6])
  abline(h=0, col="gray")
  axis(1, labels=F)
  axis(2)
  axis(2, at=0, labels="residuals", outer=T, line=1.5, tick=0, las=0)
  axis(3, at=0.5*sum(range(x1)), labels="dBS", outer=T, line=0.1, tick=0)
  box()
  
  plot(range(x1), ylims, type="n", xlab="", ylab=" ", axes=F)
  segments(x0=x1, x1=x1, y0=sumR2[,4], y1=sumR2[,8])
  segments(x0=x1, x1=x1, y0=sumR2[,5], y1=sumR2[,7], lwd=3)
  points(x1, sumR2[,6])
  abline(h=0, col="gray")
  axis(1, labels=F)
  axis(2, labels=F)
  axis(3, at=0.5*sum(range(x1)), labels="HS", line=0.1, tick=0)
  box()
  
  plot(range(x1), ylims, type="n", xlab="", ylab="", axes=F)
  segments(x0=x1, x1=x1, y0=sumR3[,4], y1=sumR3[,8])
  segments(x0=x1, x1=x1, y0=sumR3[,5], y1=sumR3[,7], lwd=3)
  points(x1, sumR3[,6])
  abline(h=0, col="gray")
  axis(1, labels=F)
  axis(2, labels=F)
  axis(3, at=0.5*sum(range(x1)), labels="SF", line=0.1, tick=0)
  box()
  
  plot(range(x1), ylims, type="n", xlab="", ylab="", axes=F)
  segments(x0=x1, x1=x1, y0=sumR4[,4], y1=sumR4[,8])
  segments(x0=x1, x1=x1, y0=sumR4[,5], y1=sumR4[,7], lwd=3)
  points(x1, sumR4[,6])
  abline(h=0, col="gray")
  axis(1, labels=F)
  axis(2, labels=F)
  axis(3, at=0.5*sum(range(x1)), labels="LM", line=0.1, tick=0)
  box()
  
#  par(mar=c(2.5, 3, 2.5, 0.25))
  #### predict v observed
  sumP1 <- summary(yhat1)  
  sumP2 <- summary(yhat2)  
  sumP3 <- summary(yhat3)  
  sumP4 <- summary(yhat4)  

  plot(x1, y1, xlab="", ylab="", axes=F)
#  segments(x0=x1, x1=x1, y0=sumP1[,4], y1=sumP1[,8], col=grey(0.4))
#  segments(x0=x1, x1=x1, y0=sumP1[,5], y1=sumP1[,7], col=grey(0.4), lwd=3)
  oo<- order(x1)
  polygon(x=c(x1[oo], rev(x1[oo])), y=c(sumP1[oo,4], rev(sumP1[oo,8])),
          col=grey(0.6))
  polygon(x=c(x1[oo], rev(x1[oo])), y=c(sumP1[oo,5], rev(sumP1[oo,7])),
          col=grey(0.4))
  lines(x1[oo], sumP1[oo,6])
  points(x1, y1)
  axis(1, labels=F)
  axis(2)
  axis(2, at=0.5*sum(range(y1)), labels=Ylab, outer=T, line=1.5, tick=0, las=0)
  box()

  plot(x1, y1, xlab="", ylab="", axes=F)
#  segments(x0=x1, x1=x1, y0=sumP2[,4], y1=sumP2[,8], col=grey(0.4))
#  segments(x0=x1, x1=x1, y0=sumP2[,5], y1=sumP2[,7], col=grey(0.4), lwd=3)
#  points(x1, sumP2[,6], col=grey(0.4), pch=16)
  polygon(x=c(x1[oo], rev(x1[oo])), y=c(sumP2[oo,4], rev(sumP2[oo,8])),
          col=grey(0.6))
  polygon(x=c(x1[oo], rev(x1[oo])), y=c(sumP2[oo,5], rev(sumP2[oo,7])),
          col=grey(0.4))
  lines(x1[oo], sumP2[oo,6])
  points(x1, y1)
  axis(1, labels=F)
  axis(2, labels=F)
  box()
  
  plot(x1, y1, xlab="", ylab="", axes=F)
#  segments(x0=x1, x1=x1, y0=sumP3[,4], y1=sumP3[,8], col=grey(0.4))
#  segments(x0=x1, x1=x1, y0=sumP3[,5], y1=sumP3[,7], col=grey(0.4), lwd=3)
#  points(x1, sumP3[,6], col=gray(0.4), pch=16)
  polygon(x=c(x1[oo], rev(x1[oo])), y=c(sumP3[oo,4], rev(sumP3[oo,8])),
          col=grey(0.6))
  polygon(x=c(x1[oo], rev(x1[oo])), y=c(sumP3[oo,5], rev(sumP3[oo,7])),
          col=grey(0.4))
  lines(x1[oo], sumP3[oo,6])
  points(x1, y1)
  axis(1, labels=F)
  axis(2, labels=F)
  box()

  plot(x1, y1, xlab="", ylab="", axes=F)
#  segments(x0=x1, x1=x1, y0=sumP3[,4], y1=sumP3[,8], col=grey(0.4))
#  segments(x0=x1, x1=x1, y0=sumP3[,5], y1=sumP3[,7], col=grey(0.4), lwd=3)
#  points(x1, sumP3[,6], col=gray(0.4), pch=16)
  polygon(x=c(x1[oo], rev(x1[oo])), y=c(sumP4[oo,4], rev(sumP4[oo,8])),
          col=grey(0.6))
  polygon(x=c(x1[oo], rev(x1[oo])), y=c(sumP4[oo,5], rev(sumP4[oo,7])),
          col=grey(0.4))
  lines(x1[oo], sumP4[oo,6])
  points(x1, y1)
  axis(1, labels=T)
  axis(2, labels=F)
  box()

#  par(mar=c(3, 3, 2, 0.25))
#  #### Bayesian p-value
#  pred1 <- rvnorm(1, mean=yhat1, sd=model1["sigma"])
#  pred2 <- rvnorm(1, mean=yhat2, sd=model2["sigma"])
#  pred3 <- rvnorm(1, mean=yhat3, sd=0.5*(model3["sigma[1]"]+model3["sigma[2]"]))
#  pred4 <- rvnorm(1, mean=yhat4, sd=model4["sigma"])
#  hist(Pr(pred1>y1), main="", xlab="Bayesian p-value")
#  hist(Pr(pred2>y1), main="", xlab="Bayesian p-value")
#  hist(Pr(pred3>y1), main="", xlab="Bayesian p-value")
#  hist(Pr(pred4>y1), main="", xlab="Bayesian p-value")
  #### change point distribution
  h1 <- hist(sims(model1["phi"]), plot=F)$density
  h2 <- hist(sims(model2["phi"]), plot=F)$density
  h3 <- hist(sims(model3["phi"]), plot=F)$density
  yrange <- range(h1, h2, h3)
  hist(sims(model1["phi"]), main="", xlab="",ylab="",
       xlim=range(x1), ylim=yrange, axes=F, col="gray", freq=F)
  axis(1, labels=T)
  axis(2, labels=T)
  axis(2, at=0.5*sum(yrange), labels="change point", outer=T,las=0,
       line=1.5, tick=F)
  box()

  hist(sims(model2["phi"]), main="", xlab="",ylab="", 
       xlim=range(x1), ylim=yrange, axes=F, col="gray", freq=F)
  axis(1, labels=T)
  axis(2, labels=F)
  box()

  hist(sims(model3["phi"]), main="", xlab="",ylab="", 
       xlim=range(x1), ylim=yrange, axes=F, col="gray", freq=F)
  axis(1, labels=T)
  axis(2, labels=F)
  box()

  plot(c(0,1),c(0,1), type="n", axes=F, xlab="", ylab="")
  text(0.5, 0.5, Txt1)
  text(0.5, 0.25, Txt2)
  mtext("log total phosphorus", side=1, outer=T, line=1.5)
invisible(list(resids1, resids2, resids3, resids4))
}

#postscript(file=paste(plotDIR, "diatom1.eps", sep="/"), height=4.5, width=6.5, horizontal=F)
diatom1.res <- threshold.plots(obsdata=diatoms[diatoms$HABITAT==subS[1] &
                  (!is.na(diatoms$GM.UTP)),], Txt1="Diatoms", Txt2=subS[1],
                y.col="Percent", x.col="GM.UTP", y.trans=logit, x.trans=log,
                simdata=simdataDiatom[[1]], Xlab="log TP", Ylab="logit % diatom")
#dev.off()

#postscript(file=paste(plotDIR, "diatom2.eps", sep="/"), height=4.5, width=6.5, horizontal=F)
diatom2.res <- threshold.plots(obsdata=diatoms[diatoms$HABITAT==subS[2] &
                  (!is.na(diatoms$GM.UTP)),], Txt1="Diatoms", Txt2=subS[2],
                y.col="Percent", x.col="GM.UTP", y.trans=logit, x.trans=log,
                simdata=simdataDiatom[[2]], Xlab="log TP", Ylab="logit % diatom")
#dev.off()

#postscript(file=paste(plotDIR, "diatom3.eps", sep="/"), height=4.5, width=6.5, horizontal=F)
diatom3.res <- threshold.plots(obsdata=diatoms[diatoms$HABITAT==subS[3] &
                  (!is.na(diatoms$GM.UTP)),], Txt1="Diatoms", Txt2=subS[3],
                y.col="Percent", x.col="GM.UTP", y.trans=logit, x.trans=log,
                simdata=simdataDiatom[[3]], Xlab="log TP", Ylab="logit % diatom")
#dev.off()

#postscript(file=paste(plotDIR, "diatom4.eps", sep="/"), height=4.5, width=6.5, horizontal=F)
diatom4.res <- threshold.plots(obsdata=diatoms[diatoms$HABITAT==subS[4] &
                  (!is.na(diatoms$GM.UTP)),], Txt1="Diatoms", Txt2=subS[4],
                y.col="Percent", x.col="GM.UTP", y.trans=logit, x.trans=log,
                simdata=simdataDiatom[[4]], Xlab="log TP", Ylab="logit % diatom")
#dev.off()

#postscript(file=paste(plotDIR, "BCD1.eps", sep="/"), height=4.5, width=6.5, horizontal=F)
bcd1.res <- threshold.plots(obsdata=macroinv[macroinv$DATE==bcdD[1]&(!is.na(macroinv$GM.UTP)),],
                y.col="BCD", x.col="GM.UTP", y.trans=function(x)x, x.trans=log,
                simdata=simdataBCD[[1]], Xlab="log TP", Ylab="BCD",
                 Txt1="BCD", Txt2=bcdD[1])
#dev.off()

#postscript(file=paste(plotDIR, "BCD2.eps", sep="/"), height=4.5, width=6.5, horizontal=F)
bcd2.res <- threshold.plots(obsdata=macroinv[macroinv$DATE==bcdD[2]&(!is.na(macroinv$GM.UTP)),],
                y.col="BCD", x.col="GM.UTP", y.trans=function(x)x, x.trans=log,
                simdata=simdataBCD[[2]], Xlab="log TP", Ylab="BCD",
                 Txt1="BCD", Txt2=bcdD[2])
#dev.off()

#postscript(file=paste(plotDIR, "BCD3.eps", sep="/"), height=4.5, width=6.5, horizontal=F)
bcd3.res <- threshold.plots(obsdata=macroinv[macroinv$DATE==bcdD[3]&(!is.na(macroinv$GM.UTP)),],
                y.col="BCD", x.col="GM.UTP", y.trans=function(x)x, x.trans=log,
                simdata=simdataBCD[[3]], Xlab="log TP", Ylab="BCD",
                 Txt1="BCD", Txt2=bcdD[3])
#dev.off()

#postscript(file=paste(plotDIR, "BCD4.eps", sep="/"), height=4.5, width=6.5, horizontal=F)
bcd4.res <- threshold.plots(obsdata=macroinv[macroinv$DATE==bcdD[4]&(!is.na(macroinv$GM.UTP)),],
                y.col="BCD", x.col="GM.UTP", y.trans=function(x)x, x.trans=log,
                simdata=simdataBCD[[4]], Xlab="log TP", Ylab="BCD",
                 Txt1="BCD", Txt2=bcdD[4])
#dev.off()

#postscript(file=paste(plotDIR, "BCD5.eps", sep="/"), height=4.5, width=6.5, horizontal=F)
bcd5.res <- threshold.plots(obsdata=macroinv[macroinv$DATE==bcdD[5]&(!is.na(macroinv$GM.UTP)),],
                y.col="BCD", x.col="GM.UTP", y.trans=function(x)x, x.trans=log,
                simdata=simdataBCD[[5]], Xlab="log TP", Ylab="BCD",
                 Txt1="BCD", Txt2=bcdD[5])
#dev.off()

#postscript(file=paste(plotDIR, "tol1.eps", sep="/"), height=4.5, width=6.5, horizontal=F)
tol1.res <- threshold.plots(obsdata=macroinv.tolerant[macroinv.tolerant$Date==tolD[1],],
                y.col="X.tolerant", x.col="TP_6mo_gm", y.trans=function(x)x,
                x.trans=log, Txt1="Tolerant Species", Txt2=tolD[1],
                simdata=simdataTOL[[1]], Xlab="log TP", Ylab="% tolerant")
#dev.off()

#postscript(file=paste(plotDIR, "tol2.eps", sep="/"), height=4.5, width=6.5, horizontal=F)
tol2.res <- threshold.plots(obsdata=macroinv.tolerant[macroinv.tolerant$Date==tolD[2],],
                y.col="X.tolerant", x.col="TP_6mo_gm", y.trans=function(x)x,
                x.trans=log, Txt1="Tolerant Species", Txt2=tolD[2],
                simdata=simdataTOL[[2]], Xlab="log TP", Ylab="% tolerant")
#dev.off()

#postscript(file=paste(plotDIR, "tol3.eps", sep="/"), height=4.5, width=6.5, horizontal=F)
tol3.res <- threshold.plots(obsdata=macroinv.tolerant[macroinv.tolerant$Date==tolD[3],],
                y.col="X.tolerant", x.col="TP_6mo_gm", y.trans=function(x)x,
                x.trans=log, Txt1="Tolerant Species", Txt2=tolD[3],
                simdata=simdataTOL[[3]], Xlab="log TP", Ylab="% tolerant")
#dev.off()

#postscript(file=paste(plotDIR, "tol4.eps", sep="/"), height=4.5, width=6.5, horizontal=F)
tol4.res <- threshold.plots(obsdata=macroinv.tolerant[macroinv.tolerant$Date==tolD[4],],
                y.col="X.tolerant", x.col="TP_6mo_gm", y.trans=function(x)x,
                x.trans=log, Txt1="Tolerant Species", Txt2=tolD[4],
                simdata=simdataTOL[[4]], Xlab="log TP", Ylab="% tolerant")
#dev.off()

#postscript(file=paste(plotDIR, "tollogit1.eps", sep="/"), height=4.5, width=6.5, horizontal=F)
tollogit1.res <- threshold.plots(obsdata=macroinv.tolerant[macroinv.tolerant$Date==tolD[1],],
                Txt1="Tolerant Species\n(logit)", Txt2=tolD[1],
                y.col="X.tolerant", x.col="TP_6mo_gm", y.trans=logit, x.trans=log,
                simdata=simdataTOL2[[1]], Xlab="log TP", Ylab="logit % tolerant")
#dev.off()

#postscript(file=paste(plotDIR, "tollogit2.eps", sep="/"), height=4.5, width=6.5,
#           horizontal=F)
tollogit2.res <- threshold.plots(obsdata=macroinv.tolerant[macroinv.tolerant$Date==tolD[2],],
                Txt1="Tolerant Species\n(logit)", Txt2=tolD[2],
                y.col="X.tolerant", x.col="TP_6mo_gm", y.trans=logit, x.trans=log,
                simdata=simdataTOL2[[2]], Xlab="log TP", Ylab="logit % tolerant")
#dev.off()

#postscript(file=paste(plotDIR, "tollogit3.eps", sep="/"), height=4.5, width=6.5,
#           horizontal=F)
tollogit3.res <- threshold.plots(obsdata=macroinv.tolerant[macroinv.tolerant$Date==tolD[3],],
                Txt1="Tolerant Species\n(logit)", Txt2=tolD[3],
                y.col="X.tolerant", x.col="TP_6mo_gm", y.trans=logit, x.trans=log,
                simdata=simdataTOL2[[3]], Xlab="log TP", Ylab="logit % tolerant")
#dev.off()
#
#postscript(file=paste(plotDIR, "tollogit4.eps", sep="/"), height=4.5, width=6.5,
#           horizontal=F)
tollogit4.res <- threshold.plots(obsdata=macroinv.tolerant[macroinv.tolerant$Date==tolD[4],],
                Txt1="Tolerant Species\n(logit)", Txt2=tolD[4],
                y.col="X.tolerant", x.col="TP_6mo_gm", y.trans=logit, x.trans=log,
                simdata=simdataTOL2[[4]], Xlab="log TP", Ylab="logit % tolerant")
#dev.off()

#postscript(file=paste(plotDIR, "senslogit1.eps", sep="/"), height=4.5, width=6.5,
#           horizontal=F)
sesnlogit1.res <- threshold.plots(obsdata=macroinv[macroinv$DATE==bcdD[1],],
                Txt1="Sensitive Species\n(logit)", Txt2=bcdD[1],
                y.col="DEP.SENS", x.col="GM.UTP", y.trans=logit, x.trans=log,
                simdata=simdataSENS[[1]], Xlab="log TP", Ylab="logit % sensitive")
#dev.off()

#postscript(file=paste(plotDIR, "senslogit2.eps", sep="/"), height=4.5, width=6.5,
#           horizontal=F)
senslogit2.res <- threshold.plots(obsdata=macroinv[macroinv$DATE==bcdD[2],],
                Txt1="Sensitive Species\n(logit)", Txt2=bcdD[2],
                y.col="DEP.SENS", x.col="GM.UTP", y.trans=logit, x.trans=log,
                simdata=simdataSENS[[2]], Xlab="log TP", Ylab="logit % sensitive")
#dev.off()

#postscript(file=paste(plotDIR, "senslogit3.eps", sep="/"), height=4.5, width=6.5,
           horizontal=F)
senslogit3.res <- threshold.plots(obsdata=macroinv[macroinv$DATE==bcdD[3],],
                Txt1="Sensitive Species\n(logit)", Txt2=bcdD[3],
                y.col="DEP.SENS", x.col="GM.UTP", y.trans=logit, x.trans=log,
                simdata=simdataSENS[[3]], Xlab="log TP", Ylab="logit % sensitive")
#dev.off()

#postscript(file=paste(plotDIR, "senslogit4.eps", sep="/"), height=4.5, width=6.5,
           horizontal=F)
senslogit4.res <- threshold.plots(obsdata=macroinv[macroinv$DATE==bcdD[4],],
                Txt1="Sensitive Species\n(logit)", Txt2=bcdD[4],
                y.col="DEP.SENS", x.col="GM.UTP", y.trans=logit, x.trans=log,
                simdata=simdataSENS[[4]], Xlab="log TP", Ylab="logit % sensitive")
#dev.off()


#postscript(file=paste(plotDIR, "sens1.eps", sep="/"), height=4.5, width=6.5, horizontal=F)
sens1.res <- threshold.plots(obsdata=macroinv[macroinv$DATE==bcdD[1],],
                Txt1="Sensitive Species", Txt2=bcdD[1],
                y.col="DEP.SENS", x.col="GM.UTP", y.trans=function(x)x, x.trans=log,
                simdata=simdataSENS2[[1]], Xlab="TP", Ylab="% sensitive")
#dev.off()

#postscript(file=paste(plotDIR, "sens2.eps", sep="/"), height=4.5, width=6.5, horizontal=F)
sens2.res <- threshold.plots(obsdata=macroinv[macroinv$DATE==bcdD[2],],
                Txt1="Sensitive Species", Txt2=bcdD[2],
                y.col="DEP.SENS", x.col="GM.UTP", y.trans=function(x)x, x.trans=log,
                simdata=simdataSENS2[[2]], Xlab="log TP", Ylab="% sensitive")
#dev.off()

#postscript(file=paste(plotDIR, "sens3.eps", sep="/"), height=4.5, width=6.5, horizontal=F)
sens3.res <- threshold.plots(obsdata=macroinv[macroinv$DATE==bcdD[3],],
                Txt1="Sensitive Species", Txt2=bcdD[3],
                y.col="DEP.SENS", x.col="GM.UTP", y.trans=function(x)x, x.trans=log,
                simdata=simdataSENS2[[3]], Xlab="log TP", Ylab="% sensitive")
#dev.off()

#postscript(file=paste(plotDIR, "sens4.eps", sep="/"), height=4.5, width=6.5, horizontal=F)
sens4.res <- threshold.plots(obsdata=macroinv[macroinv$DATE==bcdD[4],],
                Txt1="Sensitive Species", Txt2=bcdD[4],
                y.col="DEP.SENS", x.col="GM.UTP", y.trans=function(x)x, x.trans=log,
                simdata=simdataSENS2[[4]], Xlab="log TP", Ylab="% sensitive")
#dev.off()

## abrupt versus gradual changes

dbcable <- function(x, beta0, delta, phi, gamma){
  .x <- x-phi  
  return( beta0 + delta*(.x+gamma)^2*(.x > -gamma & .x <=0)/(2*gamma^2) -
    delta*(.x-gamma)^2*(.x > 0 & .x <= gamma)/(2*gamma^2) + 
      delta*(.x>0) )}
postscript(paste(plotDIR, "abrVgra.eps", sep="/"), height=2.5, width=2.5, horizontal=F)
par(mar=c(3,1,0.5,0.5), mgp=c(1.25, 0.5,0), tck= 0.01)
x <- seq(0,1,,1000)
y <- dbcable(x,beta0=0.75,delta= -0.5, phi=0.5,gamma=0.2)
x1 <- seq(0.3,0.7,,100)
y1 <- dbcable(x1,beta0=0.75,delta= -0.5,phi=0.5,gamma=0.001)

plot(x,y,type="l",xlab="", ylab="", axes=F)
lines(x1,y1,lty=2)  
lines(x[x>=0.3 & x<0.7], y[x>=0.3 & x<0.7], lwd=3, col="gray")
axis(1, at=c(0.5-0.2,0.5,0.5+0.2), labels=expression(phi-gamma, phi, phi+gamma))
box()
dev.off()

postscript(file=paste(plotDIR, "threshold3.eps", sep="/"),
           height=1.75, width=5.75, horizontal=F)
par(mfrow=c(1,4), mar=c(1.5,0.125,0.5,0.125), mgp=c(1.25,.25,0), tck=0.01)
plot(c(0,1), c(0, 1), type="n", xlab="", ylab="", axes=F)
box()
axis(1, at=0.5, label=expression(phi))
segments(x0=c(0, 0.5), x1=c(0.5,1), y0=c(0.25, 0.7), y1=c(0.5, 0.3))
text(0.05, 0.95, "dBS")
plot(c(0,1), c(0,1), type="n", xlab="", ylab="", axes=F)
box()
axis(1, at=0.5, label=expression(phi))
segments(x0=c(0, 0.5), x1=c(0.5,1), y0=c(0.5, 0.7), y1=c(0.7, 0.25))
text(0.05, 0.95, "HS")
plot(c(0,1), c(0,1), type="n", xlab="", ylab="", axes=F)
box()
axis(1, at=0.5, label=expression(phi))
segments(x0=c(0, 0.5), x1=c(0.5,1), y0=c(0.3, 0.75), y1=c(0.3, 0.75))
text(0.05, 0.95, "SF")
plot(c(0,1), c(0, 1), type="n", xlab="", ylab="", axes=F)
box()
segments(x0=0, x1=1, y0=0.3, y1=0.75)
text(0.05, 0.95, "LM")
dev.off()

#########Not used in this paper################
###############################################
## "nonparametric" deviance reduction method ##
###############################################
## the changepoint program for normal response:
chngp.nonpar <- function(infile)
{
    temp <- na.omit(infile)
    yy <- temp$y
    xx <- temp$x
    mx <- sort(unique(xx))
    m <- length(mx)
    vi <- numeric()
    vi [m] <- sum((yy - mean(yy))^2)
    for(i in 1:(m-1))
            vi[i] <- sum((yy[xx <= mx[i]] - mean(yy[xx <=
                mx[i]]))^2) + sum((yy[xx > mx[i]] - mean(
                yy[xx > mx[i]]))^2)
    chngp <- mean(mx[vi == min(vi)])
    return(chngp)
}

my.bcanon<-
function (x, nboot, theta, ..., alpha = c(0.025, 0.05, 0.1, 0.16,
    0.84, 0.9, 0.95, 0.975))
{
    n <- length(x)
    thetahat <- theta(x, ...)
    bootsam <- matrix(sample(x, size = n * nboot, replace = TRUE),
        nrow = nboot)
    thetastar <- apply(bootsam, 1, theta, ...)
    z0 <- qnorm(sum(thetastar < thetahat)/nboot)
    if (!is.finite(z0)) return(rep(NA, length(alpha)))
    u <- rep(0, n)
    for (i in 1:n) {
        u[i] <- theta(x[-i], ...)
    }
    uu <- mean(u) - u
    acc <- sum(uu * uu * uu)/(6 * (sum(uu * uu))^1.5)
    zalpha <- qnorm(alpha)
    tt <- pnorm(z0 + (z0 + zalpha)/(1 - acc * (z0 + zalpha)))
    ooo <- trunc(tt * nboot)
    if (ooo[1]>0) confpoints <- sort(thetastar)[ooo]
    else confpoints <- rep(NA, length(alpha))
    return(confpoints)
}

library(bootstrap)
n.boot <- 5000
size <- length(x1)
ylim <- NULL

thetastar1 <- bootstrap(x=1:size, nboot=n.boot, theta=function(x, infile){
    chngp.nonpar(infile[x, ])},
                       infile=data.frame(x=x1, y=y1)) $thetastar

    ci.bca1 <- my.bcanon(x=1:size, nboot=n.boot, theta=function(x, infile){
        chngp.nonpar(infile[x, ])},
                        infile=data.frame(x=x1, y=y1),
                        alpha=c(0.05, 0.95))
ylim <- c(ylim, hist(thetastar1, plot=F)$density)

thetastar2 <- bootstrap(x=1:size, nboot=n.boot, theta=function(x, infile){
    chngp.nonpar(infile[x, ])},
                       infile=data.frame(x=x2, y=y2)) $thetastar

    ci.bca2 <- my.bcanon(x=1:size, nboot=n.boot, theta=function(x, infile){
        chngp.nonpar(infile[x, ])},
                        infile=data.frame(x=x2, y=y2),
                        alpha=c(0.05, 0.95))
ylim <- c(ylim, hist(thetastar2, plot=F)$density)

thetastar3 <- bootstrap(x=1:size, nboot=n.boot, theta=function(x, infile){
    chngp.nonpar(infile[x, ])},
                       infile=data.frame(x=x3, y=y3)) $thetastar

    ci.bca3 <- my.bcanon(x=1:size, nboot=n.boot, theta=function(x, infile){
        chngp.nonpar(infile[x, ])},
                        infile=data.frame(x=x3, y=y3),
                        alpha=c(0.05, 0.95))
ylim <- c(ylim, hist(thetastar3, plot=F)$density)

thetastar4 <- bootstrap(x=1:size, nboot=n.boot, theta=function(x, infile){
    chngp.nonpar(infile[x, ])},
                       infile=data.frame(x=x4, y=y4)) $thetastar

    ci.bca4 <- my.bcanon(x=1:size, nboot=n.boot, theta=function(x, infile){
        chngp.nonpar(infile[x, ])},
                        infile=data.frame(x=x4, y=y4),
                        alpha=c(0.05, 0.95))
ylim <- c(ylim, hist(thetastar4, plot=F)$density)

par(mfrow=c(2,2), mar=c(3,3,3,0.25), mgp=c(1.5,0.5,0), tck= -0.01)
hist(thetastar1, col="gray", border=F, xlim=c(0,1), main="general",
     xlab="x", axes=F, ylab="", ylim=c(0, max(ylim)), prob=T)
segments(x0=quantile(thetastar1, prob=c(0.025, 0.975)),
         x1=quantile(thetastar1, prob=c(0.025, 0.975)),
         y0=c(-0.05,-0.050)*max(ylim), y1=rep(0.025*max(ylim), 2), lwd=2)
axis(1, at=0.5)
axis(1)

hist(thetastar2, col="gray", border=F, xlim=c(0,1), main="hockey stick",
     xlab="x", axes=F, ylab="", ylim=c(0, max(ylim)), prob=T)
segments(x0=quantile(thetastar2, prob=c(0.025, 0.975)),
         x1=quantile(thetastar2, prob=c(0.025, 0.975)),
         y0=c(-0.05,-0.050)*max(ylim), y1=rep(0.025*max(ylim), 2), lwd=2)
axis(1, at=0.5)
axis(1)

hist(thetastar3, col="gray", border=F, xlim=c(0,1), main="step function",
     xlab="x", axes=F, ylab="", ylim=c(0, max(ylim)), prob=T)
segments(x0=quantile(thetastar3, prob=c(0.025, 0.975)),
         x1=quantile(thetastar3, prob=c(0.025, 0.975)),
         y0=c(-0.05,-0.050)*max(ylim), y1=rep(0.025*max(ylim), 2), lwd=2)
axis(1, at=0.5)
axis(1)

hist(thetastar4, col="gray", border=F, xlim=c(0,1), main="linear",
     xlab="x", axes=F, ylab="", ylim=c(0, max(ylim)), prob=T)
segments(x0=quantile(thetastar4, prob=c(0.025, 0.975)),
         x1=quantile(thetastar4, prob=c(0.025, 0.975)),
         y0=c(-0.05,-0.050)*max(ylim), y1=rep(0.025*max(ylim), 2), lwd=2)
axis(1, at=0.5)
axis(1)

mu.ci <-
rbind(temp[1:3,],
c(chngp.nonpar(data.frame(x=x1, y=y1)),quantile(thetastar1, prob=c(0.025, 0.975))),
      temp[4:6,],
c(chngp.nonpar(data.frame(x=x2, y=y2)),
quantile(thetastar2, prob=c(0.025, 0.975))),
      temp[7:9,],
c(chngp.nonpar(data.frame(x=x3, y=y3)),
quantile(thetastar3, prob=c(0.025, 0.975))),
      temp[10:12,],
      c(chngp.nonpar(data.frame(x=x4, y=y4)),
quantile(thetastar4, prob=c(0.025, 0.975))))

write(round(t(mu.ci[1:4,]), 3),
      file="thethaG.txt", sep="&", ncolumns=dim(mu.ci)[2])
write(round(t(mu.ci[5:8,]), 3),
      file="thethaH.txt", sep="&", ncolumns=dim(mu.ci)[2])
write(round(t(mu.ci[9:12,]), 3),
      file="thethaS.txt", sep="&", ncolumns=dim(mu.ci)[2])
write(round(t(mu.ci[13:16,]), 3),
      file="thethaL.txt", sep="&", ncolumns=dim(mu.ci)[2])

## residual plot

x <- runif(40)
y <- 2+3*x+rnorm(40, 0, 0.25)
m <- c(mean(y[x<0.5]), mean(y[x>0.5]))
xx <- c(x,0.5001, 0.4999)
oo <- order(xx)
resSF <- y-m[as.numeric(x>0.5)+1]
postscript(paste(plotDIR, "residSF.eps",sep="/"),
           height=3.5,width=2.5,horizontal=F)
par(mfrow=c(3,1), oma=c(3,3,1,3), mar=c(0.125,0.125,0.125,0.125),
    mgp=c(1.25,0.25,0), las=1, tck=0.01)
plot(x,y, xlab="", ylab="", axes=F)
box()
axis(2)
#abline(v=0.5)
abline(lsfit(x, y), lty=2, col="gray")
lines(xx[oo], m[as.numeric(xx>0.5)+1][oo], col="gray")
plot(resid(lm(y~x))~x, xlab="",ylab="", axes=F)
box()
axis(4)
plot(x, resSF, xlab="x", ylab="", axes=F)
box()
axis(1)
axis(2)
dev.off()

## residual Q-Q norm plots:
qqnorm.rv <- function(x.rv, Prob=c(0.01,0.025,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.975,0.99), ...){
    plot.rv(qnorm(Prob), quantile(x.rv, prob=Prob), ...)
    invisible()
}

residual.qqplots <-function(obsdata=diatoms[diatoms$HABITAT==subS[2] & (!is.na(diatoms$GM.UTP)),],
                            Txt1="Diatoms", Txt2=subS[2],
                            y.col="Percent", x.col="GM.UTP",
                            y.trans=logit, x.trans=log,
                            simdata=simdataDiatom[[2]],
                            Xlab="log TP", Ylab="logit % diatom"){
    ## extracting model results
    model1 <- rvsims(simdata[[1]]) ## gen
    model2 <- rvsims(simdata[[2]]) ## hockey
    model3 <- rvsims(simdata[[3]]) ## step
    model4 <- rvsims(simdata[[4]]) ## linear
    print(1)
    ## extracting data
    y1 <- y.trans(obsdata[,y.col])
    x1 <- x.trans(obsdata[,x.col])
    resids1<-y1-(yhat1<-(model1["beta[1]"]+ model1["delta[1]"]*(x1>=model1["phi"]))+ (model1["beta[2]"]+model1["delta[2]"]*(x1>=model1["phi"]))*(x1-model1["phi"]))
    resids.range <- range(summary(quantile.rv(resids1, prob=c(0.01,0.99)))[c(1,2), c(4, 10)])
    print(2)
    resids2<-y1-(yhat2<-model2["beta[1]"] + (model2["beta[2]"]+model2["delta"]*(x1>=model2["phi"]))*(x1-model2["phi"]))
    print(3)
    resids.range <- c(resids.range, range(summary(quantile.rv(resids2, prob=c(0.01,0.99)))[c(1,2), c(4, 10)]))
    resids3<-y1-(yhat3<-model3["beta"] + model3["delta"]*(x1>=model3["phi"]))
    print(4)
    resids.range <- c(resids.range, range(summary(quantile.rv(resids3, prob=c(0.01,0.99)))[c(1,2), c(4, 10)]))
    resids4<-y1-(yhat4<-model4["beta[1]"] + model4["beta[2]"]*x1)
    print(5)
    ## plots
    resids.range <- c(resids.range, range(summary(quantile.rv(resids4, prob=c(0.025,0.975)))[c(1,2), c(4, 10)]))
    if(range(resids.range)[2]>1){
        Yloc <- floor(range(resids.range)[2])
    } else {
        Yloc <- range(resids.range)[2]
    }
    Yloc2 <- range(resids.range)[1]+0.2*diff(range(resids.range))
    par(mfrow=c(1,4), oma=c(3, 3, 2, 0.25),
        mar=c(3, 0.125, 3, 0.125),
        mgp=c(1.25,0.125,0.), las=1, tck=0.01)
    qqnorm.rv(resids1, ylab="residual quantiles",xlab="",
              axes=F, ylim=range(resids.range))
    axis(1)
    axis(2)
    box()
    text(-1,Yloc, "dBS Model")
    qqnorm.rv(resids2, ylab="residual quantiles",xlab="",
              axes=F, ylim=range(resids.range))
    axis(3)
    box()
    text(-1,Yloc, "BS Model")
    qqnorm.rv(resids3, ylab="residual quantiles",xlab="",
              axes=F, ylim=range(resids.range))
    axis(1)
    box()
    text(-1,Yloc, "SF Model")
    qqnorm.rv(resids4, ylab="residual quantiles",xlab="",
              axes=F, ylim=range(resids.range))
    axis(3)
    axis(4)
    box()
    text(-1,Yloc, "Linear Model")
    text(1, Yloc2, paste(Txt1, Txt2, sep="\n"))
    invisible()
}

postscript(file=paste(plotDIR, "diatom1resQQ.eps", sep="/"), height=3.1, width=6.5, horizontal=F)
residual.qqplots(obsdata=diatoms[diatoms$HABITAT==subS[1] &
                  (!is.na(diatoms$GM.UTP)),], Txt1="% Diatoms", Txt2=subS[1],
                y.col="Percent", x.col="GM.UTP", y.trans=logit, x.trans=log,
                simdata=simdataDiatom[[1]], Xlab="log TP", Ylab="logit % diatom")
dev.off()

postscript(file=paste(plotDIR, "diatom2resQQ.eps", sep="/"), height=3.1, width=6.5, horizontal=F)
residual.qqplots(obsdata=diatoms[diatoms$HABITAT==subS[2] &
                  (!is.na(diatoms$GM.UTP)),], Txt1="% Diatoms", Txt2=subS[2],
                y.col="Percent", x.col="GM.UTP", y.trans=logit, x.trans=log,
                simdata=simdataDiatom[[2]], Xlab="log TP", Ylab="logit % diatom")
dev.off()

postscript(file=paste(plotDIR, "diatom3resQQ.eps", sep="/"), height=3.1, width=6.5, horizontal=F)
residual.qqplots(obsdata=diatoms[diatoms$HABITAT==subS[3] &
                  (!is.na(diatoms$GM.UTP)),], Txt1="% Diatoms", Txt2=subS[3],
                y.col="Percent", x.col="GM.UTP", y.trans=logit, x.trans=log,
                simdata=simdataDiatom[[3]], Xlab="log TP", Ylab="logit % diatom")
dev.off()

postscript(file=paste(plotDIR, "diatom4resQQ.eps", sep="/"), height=3.1, width=6.5, horizontal=F)
residual.qqplots(obsdata=diatoms[diatoms$HABITAT==subS[4] &
                  (!is.na(diatoms$GM.UTP)),], Txt1="% Diatoms", Txt2=subS[4],
                y.col="Percent", x.col="GM.UTP", y.trans=logit, x.trans=log,
                simdata=simdataDiatom[[4]], Xlab="log TP", Ylab="logit % diatom")
dev.off()

postscript(file=paste(plotDIR, "BCD1resQQ.eps", sep="/"), height=3.1, width=6.5, horizontal=F)
residual.qqplots(obsdata=macroinv[macroinv$DATE==bcdD[1]&(!is.na(macroinv$GM.UTP)),],
                y.col="BCD", x.col="GM.UTP", y.trans=function(x)x, x.trans=log,
                simdata=simdataBCD[[1]], Xlab="log TP", Ylab="BCD",
                 Txt1="BCD", Txt2=bcdD[1])
dev.off()

postscript(file=paste(plotDIR, "BCD2resQQ.eps", sep="/"), height=3.1, width=6.5, horizontal=F)
residual.qqplots(obsdata=macroinv[macroinv$DATE==bcdD[2]&(!is.na(macroinv$GM.UTP)),],
                y.col="BCD", x.col="GM.UTP", y.trans=function(x)x, x.trans=log,
                simdata=simdataBCD[[2]], Xlab="log TP", Ylab="BCD",
                 Txt1="BCD", Txt2=bcdD[2])
dev.off()

postscript(file=paste(plotDIR, "BCD3resQQ.eps", sep="/"), height=3.1, width=6.5, horizontal=F)
residual.qqplots(obsdata=macroinv[macroinv$DATE==bcdD[3]&(!is.na(macroinv$GM.UTP)),],
                y.col="BCD", x.col="GM.UTP", y.trans=function(x)x, x.trans=log,
                simdata=simdataBCD[[3]], Xlab="log TP", Ylab="BCD",
                 Txt1="BCD", Txt2=bcdD[3])
dev.off()

postscript(file=paste(plotDIR, "BCD4resQQ.eps", sep="/"), height=3.1, width=6.5, horizontal=F)
residual.qqplots(obsdata=macroinv[macroinv$DATE==bcdD[4]&(!is.na(macroinv$GM.UTP)),],
                y.col="BCD", x.col="GM.UTP", y.trans=function(x)x, x.trans=log,
                simdata=simdataBCD[[4]], Xlab="log TP", Ylab="BCD",
                 Txt1="BCD", Txt2=bcdD[4])
dev.off()

postscript(file=paste(plotDIR, "BCD5resQQ.eps", sep="/"), height=3.1, width=6.5, horizontal=F)
residual.qqplots(obsdata=macroinv[macroinv$DATE==bcdD[5]&(!is.na(macroinv$GM.UTP)),],
                y.col="BCD", x.col="GM.UTP", y.trans=function(x)x, x.trans=log,
                simdata=simdataBCD[[5]], Xlab="log TP", Ylab="BCD",
                 Txt1="BCD", Txt2=bcdD[5])
dev.off()

postscript(file=paste(plotDIR, "tol1resQQ.eps", sep="/"), height=3.1, width=6.5, horizontal=F)
residual.qqplots(obsdata=macroinv.tolerant[macroinv.tolerant$Date==tolD[1],],
                y.col="X.tolerant", x.col="TP_6mo_gm", y.trans=function(x)x,
                x.trans=log, Txt1="% Tolerant", Txt2=tolD[1],
                simdata=simdataTOL[[1]], Xlab="log TP", Ylab="% tolerant")
dev.off()

postscript(file=paste(plotDIR, "tol2resQQ.eps", sep="/"), height=3.1, width=6.5, horizontal=F)
residual.qqplots(obsdata=macroinv.tolerant[macroinv.tolerant$Date==tolD[2],],
                y.col="X.tolerant", x.col="TP_6mo_gm", y.trans=function(x)x,
                x.trans=log, Txt1="% Tolerant", Txt2=tolD[2],
                simdata=simdataTOL[[2]], Xlab="log TP", Ylab="% tolerant")
dev.off()

postscript(file=paste(plotDIR, "tol3resQQ.eps", sep="/"), height=3.1, width=6.5, horizontal=F)
residual.qqplots(obsdata=macroinv.tolerant[macroinv.tolerant$Date==tolD[3],],
                y.col="X.tolerant", x.col="TP_6mo_gm", y.trans=function(x)x,
                x.trans=log, Txt1="% Tolerant", Txt2=tolD[3],
                simdata=simdataTOL[[3]], Xlab="log TP", Ylab="% tolerant")
dev.off()

postscript(file=paste(plotDIR, "tol4resQQ.eps", sep="/"), height=3.1, width=6.5, horizontal=F)
residual.qqplots(obsdata=macroinv.tolerant[macroinv.tolerant$Date==tolD[4],],
                y.col="X.tolerant", x.col="TP_6mo_gm", y.trans=function(x)x,
                x.trans=log, Txt1="% Tolerant", Txt2=tolD[4],
                simdata=simdataTOL[[4]], Xlab="log TP", Ylab="% tolerant")
dev.off()

postscript(file=paste(plotDIR, "tollogit1resQQ.eps", sep="/"), height=3.1, width=6.5, horizontal=F)
residual.qqplots(obsdata=macroinv.tolerant[macroinv.tolerant$Date==tolD[1],],
                Txt1="logit % Tolerant", Txt2=tolD[1],
                y.col="X.tolerant", x.col="TP_6mo_gm", y.trans=logit, x.trans=log,
                simdata=simdataTOL2[[1]], Xlab="log TP", Ylab="logit % tolerant")
dev.off()

postscript(file=paste(plotDIR, "tollogit2resQQ.eps", sep="/"), height=3.1, width=6.5,
           horizontal=F)
residual.qqplots(obsdata=macroinv.tolerant[macroinv.tolerant$Date==tolD[2],],
                Txt1="logit % Tolerant", Txt2=tolD[2],
                y.col="X.tolerant", x.col="TP_6mo_gm", y.trans=logit, x.trans=log,
                simdata=simdataTOL2[[2]], Xlab="log TP", Ylab="logit % tolerant")
dev.off()

postscript(file=paste(plotDIR, "tollogit3resQQ.eps", sep="/"), height=3.1, width=6.5,
           horizontal=F)
residual.qqplots(obsdata=macroinv.tolerant[macroinv.tolerant$Date==tolD[3],],
                Txt1="logit % Tolerant", Txt2=tolD[3],
                y.col="X.tolerant", x.col="TP_6mo_gm", y.trans=logit, x.trans=log,
                simdata=simdataTOL2[[3]], Xlab="log TP", Ylab="logit % tolerant")
dev.off()

postscript(file=paste(plotDIR, "tollogit4resQQ.eps", sep="/"), height=3.1, width=6.5,
           horizontal=F)
residual.qqplots(obsdata=macroinv.tolerant[macroinv.tolerant$Date==tolD[4],],
                Txt1="logit % Tolerant", Txt2=tolD[4],
                y.col="X.tolerant", x.col="TP_6mo_gm", y.trans=logit, x.trans=log,
                simdata=simdataTOL2[[4]], Xlab="log TP", Ylab="logit % tolerant")
dev.off()

postscript(file=paste(plotDIR, "senslogit1resQQ.eps", sep="/"), height=3.1, width=6.5,
           horizontal=F)
residual.qqplots(obsdata=macroinv[macroinv$DATE==bcdD[1],],
                Txt1="logit % Sensitive", Txt2=bcdD[1],
                y.col="DEP.SENS", x.col="GM.UTP", y.trans=logit, x.trans=log,
                simdata=simdataSENS[[1]], Xlab="log TP", Ylab="logit % sensitive")
dev.off()

postscript(file=paste(plotDIR, "senslogit2resQQ.eps", sep="/"), height=3.1, width=6.5,
           horizontal=F)
residual.qqplots(obsdata=macroinv[macroinv$DATE==bcdD[2],],
                Txt1="logit % Sensitive", Txt2=bcdD[2],
                y.col="DEP.SENS", x.col="GM.UTP", y.trans=logit, x.trans=log,
                simdata=simdataSENS[[2]], Xlab="log TP", Ylab="logit % sensitive")
dev.off()

postscript(file=paste(plotDIR, "senslogit3resQQ.eps", sep="/"), height=3.1, width=6.5,
           horizontal=F)
residual.qqplots(obsdata=macroinv[macroinv$DATE==bcdD[3],],
                Txt1="logit % Sensitive", Txt2=bcdD[3],
                y.col="DEP.SENS", x.col="GM.UTP", y.trans=logit, x.trans=log,
                simdata=simdataSENS[[3]], Xlab="log TP", Ylab="logit % sensitive")
dev.off()

postscript(file=paste(plotDIR, "senslogit4resQQ.eps", sep="/"), height=3.1, width=6.5,
           horizontal=F)
residual.qqplots(obsdata=macroinv[macroinv$DATE==bcdD[4],],
                Txt1="logit % Sensitive", Txt2=bcdD[4],
                y.col="DEP.SENS", x.col="GM.UTP", y.trans=logit, x.trans=log,
                simdata=simdataSENS[[4]], Xlab="log TP", Ylab="logit % sensitive")
dev.off()


postscript(file=paste(plotDIR, "sens1resQQ.eps", sep="/"), height=3.1, width=6.5, horizontal=F)
residual.qqplots(obsdata=macroinv[macroinv$DATE==bcdD[1],],
                Txt1="% Sensitive", Txt2=bcdD[1],
                y.col="DEP.SENS", x.col="GM.UTP", y.trans=function(x)x, x.trans=log,
                simdata=simdataSENS2[[1]], Xlab="TP", Ylab="% sensitive")
dev.off()

postscript(file=paste(plotDIR, "sens2resQQ.eps", sep="/"), height=3.1, width=6.5, horizontal=F)
residual.qqplots(obsdata=macroinv[macroinv$DATE==bcdD[2],],
                Txt1="% Sensitive", Txt2=bcdD[2],
                y.col="DEP.SENS", x.col="GM.UTP", y.trans=function(x)x, x.trans=log,
                simdata=simdataSENS2[[2]], Xlab="log TP", Ylab="% sensitive")
dev.off()

postscript(file=paste(plotDIR, "sens3resQQ.eps", sep="/"), height=3.1, width=6.5, horizontal=F)
residual.qqplots(obsdata=macroinv[macroinv$DATE==bcdD[3],],
                Txt1="% Sensitive", Txt2=bcdD[3],
                y.col="DEP.SENS", x.col="GM.UTP", y.trans=function(x)x, x.trans=log,
                simdata=simdataSENS2[[3]], Xlab="log TP", Ylab="% sensitive")
dev.off()

postscript(file=paste(plotDIR, "sens4resQQ.eps", sep="/"), height=3.1, width=6.5, horizontal=F)
residual.qqplots(obsdata=macroinv[macroinv$DATE==bcdD[4],],
                Txt1="% Sensitive", Txt2=bcdD[4],
                y.col="DEP.SENS", x.col="GM.UTP", y.trans=function(x)x, x.trans=log,
                simdata=simdataSENS2[[4]], Xlab="log TP", Ylab="% sensitive")
dev.off()

