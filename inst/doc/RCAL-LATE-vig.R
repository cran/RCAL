### R code from vignette source 'RCAL-LATE-vig.Rnw'

###################################################
### code chunk number 1: read.iv.data
###################################################
library(RCAL)

data(simu.iv.data)


###################################################
### code chunk number 2: RCAL-LATE-vig.Rnw:80-81
###################################################
simu.iv.data[1:10, 1:6]


###################################################
### code chunk number 3: RCAL-LATE-vig.Rnw:87-95
###################################################
n <- dim(simu.iv.data)[1]
p <- 100 # include the first 100 covariates due to CRAN time constraint

y <- simu.iv.data[,1]
tr <- simu.iv.data[,2]
iv <- simu.iv.data[,3]
x <- simu.iv.data[,3+1:p]
x <- scale(x)


###################################################
### code chunk number 4: RCAL-LATE-vig.Rnw:104-109
###################################################
par(mfrow=c(3,2))
par(mar=c(4,4,2,2))
for (j in 1:6) {
  boxplot(x[,j] ~ iv, ylab=paste("covariate x", j, sep=""), xlab="instrument")
}


###################################################
### code chunk number 5: RCAL-LATE-vig.Rnw:118-123
###################################################
par(mfrow=c(3,2))
par(mar=c(4,4,2,2))
for (j in 1:6) {
  boxplot(x[,j] ~ tr, ylab=paste("covariate x", j, sep=""), xlab="treatment")
}


###################################################
### code chunk number 6: RCAL-LATE-vig.Rnw:132-137
###################################################
par(mfrow=c(3,2))
par(mar=c(4,4,2,2))
for (j in 1:6) {
  plot(x[tr==1,j], y[tr==1], ylab="y", xlab=paste("covariate x", j, sep=""))
}


###################################################
### code chunk number 7: RCAL-LATE-vig.Rnw:144-152
###################################################
pi    <- mean(iv)
del   <- iv/pi-(1-iv)/(1-pi);
del.y <- mean(del*y);
del.tr<- mean(del*tr);
(w<- del.y/del.tr)  # point estimate

g  <- mean((-iv/pi^2-(1-iv)/(1-pi)^2)*(y-w*tr))
sqrt(mean((del*(y-w*tr)-g*(pi-iv))^2)/del.tr^2/n) # standard error


###################################################
### code chunk number 8: RCAL-LATE-vig.Rnw:193-214
###################################################
## regularized calibrated estimation
RNGversion('3.5.0')
set.seed(0) # this affects random split of data in cross validation
late.cv.rcal <- 
late.regu.cv(fold=5*c(1,1,1), nrho=(1+10)*c(1,1,1), rho.seq=NULL,
y, tr, iv, fx=x, gx=x, hx=x, arm=2, d1=1, d2=3, ploss="cal", yloss="gaus")

matrix(unlist(late.cv.rcal$est), ncol=2, byrow=TRUE,
dimnames=list(c("ipw", "or", "est", "var", "ze",
"late.est", "late.var", "late.ze"), c("theta1", "theta0")))

## For comparison, we use the same set of models 
## in regularized maximum likelihood estimation
set.seed(0)
late.cv.rml <- 
late.regu.cv(fold=5*c(1,1,1), nrho=(1+10)*c(1,1,1), rho.seq=NULL,
y, tr, iv, fx=x, gx=x, hx=x, arm=2, d1=1, d2=3, ploss="ml", yloss="gaus")

matrix(unlist(late.cv.rml$est), ncol=2, byrow=TRUE,
dimnames=list(c("ipw", "or", "est", "var", "ze",
"late.est", "late.var", "late.ze"), c("theta1", "theta0")))


###################################################
### code chunk number 9: RCAL-LATE-vig.Rnw:231-305
###################################################
## Computation of standardized calibration differences 
##               within iv=1 group
fips.raw <- rep(mean(iv), n)   #constant propensity scores
fips1.cv.rcal<-late.cv.rcal$mfp[,2]
fips1.cv.rml <-late.cv.rml$mfp[,2]

check1.raw <- mn.ipw(x, iv, fips.raw)
check1.cv.rcal <- mn.ipw(x, iv, fips1.cv.rcal)
check1.cv.rml <- mn.ipw(x, iv, fips1.cv.rml)

par(mfrow=c(2,2))
par(mar=c(4,4,2,2))
plot(check1.raw$est, xlim=c(0,p), ylim=c(-.3,.3), 
     xlab="", ylab="std diff", main="Raw") 
abline(h=0)

plot(check1.cv.rml$est, xlim=c(0,p), ylim=c(-.3,.3), 
     xlab="", ylab="std diff", main="RML") 
abline(h=0)
abline(h=max(abs(check1.cv.rml$est)) *c(-1,1), lty=2)

nz1.rml <- which(late.cv.rml$ips[[2]]$sel.bet[,1]!=0)-1
text(nz1.rml, -0.3, "x", cex=1)
length(nz1.rml)

plot(check1.cv.rcal$est, xlim=c(0,p), ylim=c(-.3,.3), 
     xlab="", ylab="std diff", main="RCAL") 
abline(h=0)
abline(h=max(abs(check1.cv.rcal$est)) *c(-1,1), lty=2)

nz1.rcal <- which(late.cv.rcal$ips[[2]]$sel.bet[,1]!=0)-1
text(nz1.rcal, -0.3, "x", cex=1)
length(nz1.rcal)

plot(fips1.cv.rml[iv==1], fips1.cv.rcal[iv==1], xlim=c(0,1), ylim=c(0,1), 
     xlab="RML", ylab="RCAL", main="fitted probabilities")
abline(c(0,1))

## Computation of standardized calibration differences 
##            within iv=0 group
fips0.cv.rcal<-late.cv.rcal$mfp[,1]
fips0.cv.rml <-late.cv.rml$mfp[,1]

check0.raw <- mn.ipw(x, 1-iv, fips.raw)
check0.cv.rcal <- mn.ipw(x, 1-iv, fips0.cv.rcal)
check0.cv.rml <- mn.ipw(x, 1-iv, fips0.cv.rml)

par(mfrow=c(2,2))
par(mar=c(4,4,2,2))
plot(check0.raw$est, xlim=c(0,p), ylim=c(-.3,.3), 
     xlab="", ylab="std diff", main="Raw") 
abline(h=0)

plot(check0.cv.rml$est, xlim=c(0,p), ylim=c(-.3,.3), 
     xlab="", ylab="std diff", main="RML") 
abline(h=0)
abline(h=max(abs(check0.cv.rml$est)) *c(-1,1), lty=2)

nz0.rml <- which(late.cv.rml$ips[[2]]$sel.bet[,1]!=0)-1
text(nz0.rml, -0.3, "x", cex=1)
length(nz0.rml)

plot(check0.cv.rcal$est, xlim=c(0,p), ylim=c(-.3,.3), 
     xlab="", ylab="std diff", main="RCAL") 
abline(h=0)
abline(h=max(abs(check0.cv.rcal$est)) *c(-1,1), lty=2)

nz0.rcal <- which(late.cv.rcal$ips[[1]]$sel.bet[,1]!=0)-1
text(nz0.rcal, -0.3, "x", cex=1)
length(nz0.rcal)

plot(fips0.cv.rml[iv==0], fips0.cv.rcal[iv==0], xlim=c(0,1), ylim=c(0,1), 
     xlab="RML", ylab="RCAL", main="fitted probabilities")
abline(c(0,1))


###################################################
### code chunk number 10: RCAL-LATE-vig.Rnw:312-338
###################################################

par(mfrow=c(2,2))
par(mar=c(4,4,2,2))
plot(check1.raw$est, xlim=c(0,p), ylim=c(-.3,.3), 
     xlab="", ylab="std diff", main="Raw") 
abline(h=0)

plot(check1.cv.rml$est, xlim=c(0,p), ylim=c(-.3,.3), 
     xlab="", ylab="std diff", main="RML") 
abline(h=0)
abline(h=max(abs(check1.cv.rml$est)) *c(-1,1), lty=2)

nz=which(late.cv.rml$ips[[2]]$sel.bet[,1]!=0)-1
text(nz, -0.3, "x", cex=1)

plot(check1.cv.rcal$est, xlim=c(0,p), ylim=c(-.3,.3), 
     xlab="", ylab="std diff", main="RCAL") 
abline(h=0)
abline(h=max(abs(check1.cv.rcal$est)) *c(-1,1), lty=2)

nz=which(late.cv.rcal$ips[[2]]$sel.bet[,1]!=0)-1
text(nz, -0.3, "x", cex=1)

plot(fips1.cv.rml[iv==1], fips1.cv.rcal[iv==1], xlim=c(0,1), ylim=c(0,1), 
     xlab="RML", ylab="RCAL", main="fitted probabilities")
abline(c(0,1))


###################################################
### code chunk number 11: RCAL-LATE-vig.Rnw:347-372
###################################################
par(mfrow=c(2,2))
par(mar=c(4,4,2,2))
plot(check0.raw$est, xlim=c(0,p), ylim=c(-.3,.3), 
     xlab="", ylab="std diff", main="Raw") 
abline(h=0)

plot(check0.cv.rml$est, xlim=c(0,p), ylim=c(-.3,.3), 
     xlab="", ylab="std diff", main="RML") 
abline(h=0)
abline(h=max(abs(check0.cv.rml$est)) *c(-1,1), lty=2)

nz=which(late.cv.rml$ips[[2]]$sel.bet[,1]!=0)-1
text(nz, -0.3, "x", cex=1)

plot(check0.cv.rcal$est, xlim=c(0,p), ylim=c(-.3,.3), 
     xlab="", ylab="std diff", main="RCAL") 
abline(h=0)
abline(h=max(abs(check0.cv.rcal$est)) *c(-1,1), lty=2)

nz=which(late.cv.rcal$ips[[1]]$sel.bet[,1]!=0)-1
text(nz, -0.3, "x", cex=1)

plot(fips0.cv.rml[iv==0], fips0.cv.rcal[iv==0], xlim=c(0,1), ylim=c(0,1), 
     xlab="RML", ylab="RCAL", main="fitted probabilities")
abline(c(0,1))


###################################################
### code chunk number 12: RCAL-LATE-vig.Rnw:380-434
###################################################
set.seed(0)
late.path.rcal <- 
late.regu.path(fold=5*c(1,1,1), nrho=(1+10)*c(1,1,1), rho.seq=NULL,
y, tr, iv, fx=x, gx=x, hx=x, arm=1, d1=1, d2=3, ploss="cal", yloss="gaus")

set.seed(0)
late.path.rml <- 
late.regu.path(fold=5*c(1,1,1), nrho=(1+10)*c(1,1,1), rho.seq=NULL,
y, tr, iv, fx=x, gx=x, hx=x, arm=1, d1=1, d2=3, ploss="ml", yloss="gaus")
## Computation of standardized calibration differences 
##    along regularization path within iv=1 group
fips1.path.rcal <- late.path.rcal$mfp[[2]]
mdiff1.path.rcal <- rep(NA, dim(fips1.path.rcal)[2])
rvar1.path.rcal <- rep(NA, dim(fips1.path.rcal)[2])
for (j in 1:dim(fips1.path.rcal)[2]) {
  check1.path.rcal <- mn.ipw(x, iv, fips1.path.rcal[,j])
  mdiff1.path.rcal[j] <- max(abs(check1.path.rcal$est))
  rvar1.path.rcal[j] <- 
  var(1/fips1.path.rcal[iv==1,j])/mean(1/fips1.path.rcal[iv==1,j])^2
}

fips1.path.rml <- late.path.rml$mfp[[2]]
mdiff1.path.rml <- rep(NA, dim(fips1.path.rml)[2])
rvar1.path.rml <- rep(NA, dim(fips1.path.rml)[2])
for (j in 1:dim(fips1.path.rml)[2]) {
  check1.path.rml <- mn.ipw(x, iv, fips1.path.rml[,j])
  mdiff1.path.rml[j] <- max(abs(check1.path.rml$est))
  rvar1.path.rml[j] <- 
  var(1/fips1.path.rml[iv==1,j])/mean(1/fips1.path.rml[iv==1,j])^2
}

par(mfrow=c(1,2))
par(mar=c(4,4,2,2))
plot(late.path.rml$ips[[2]]$nz.all, mdiff1.path.rml, 
     xlim=c(0,p), ylim=c(0,.4), xlab="# nonzero", ylab="std diff")
lines(late.path.rml$ips[[2]]$nz.all[!late.path.rml$ips[[2]]$non.conv],
       mdiff1.path.rml, lty=3)

points(late.path.rcal$ips[[2]]$nz.all[!late.path.rcal$ips[[2]]$non.conv],
       mdiff1.path.rcal, pch=4)
lines(late.path.rcal$ips[[2]]$nz.all[!late.path.rcal$ips[[2]]$non.conv],
       mdiff1.path.rcal, lty=3)

legend(120,.4, c("RML","RCAL"), pch=c(1,4), cex=.6)

#
plot(rvar1.path.rml, mdiff1.path.rml, 
     xlim=c(0,1), ylim=c(0,.4), xlab="rel var", ylab="std diff")
lines(rvar1.path.rml, mdiff1.path.rml, lty=3)

points(rvar1.path.rcal, mdiff1.path.rcal, pch=4)
lines(rvar1.path.rcal, mdiff1.path.rcal, lty=3)

legend(0.6,.4, c("RML","RCAL"), pch=c(1,4), cex=.6)


###################################################
### code chunk number 13: RCAL-LATE-vig.Rnw:437-483 (eval = FALSE)
###################################################
## set.seed(0)
## late.path.rcal <- 
## late.regu.path(fold=5*c(1,1,1), nrho=(1+10)*c(1,1,1), rho.seq=NULL,
## y, tr, iv, fx=x, gx=x, hx=x, arm=0, d1=1, d2=3, ploss="cal", yloss="gaus")
## 
## set.seed(0)
## late.path.rml <- 
## late.regu.path(fold=5*c(1,1,1), nrho=(1+10)*c(1,1,1), rho.seq=NULL,
## y, tr, iv, fx=x, gx=x, hx=x, arm=0, d1=1, d2=3, ploss="ml", yloss="gaus")
## ## Computation of standardized calibration differences 
## ##    along regularization path within iv=0 group
## fips0.path.rcal <- late.path.rcal$mfp[[1]]
## mdiff0.path.rcal <- rep(NA, dim(fips0.path.rcal)[2])
## rvar0.path.rcal <- rep(NA, dim(fips0.path.rcal)[2])
## for (j in 1:dim(fips0.path.rcal)[2]) {
##   check0.path.rcal <- mn.ipw(x, 1-iv, fips0.path.rcal[,j])
##   mdiff0.path.rcal[j] <- max(abs(check0.path.rcal$est))
##   rvar0.path.rcal[j] <- 
##   var(1/fips0.path.rcal[iv==0,j])/mean(1/fips0.path.rcal[iv==0,j])^2
## }
## fips0.path.rml <- late.path.rml$mfp[[1]]
## mdiff0.path.rml <- rep(NA, dim(fips0.path.rml)[2])
## rvar0.path.rml <- rep(NA, dim(fips0.path.rml)[2])
## for (j in 1:dim(fips0.path.rml)[2]) {
##   check0.path.rml <- mn.ipw(x, 1-iv, fips0.path.rml[,j])
##   mdiff0.path.rml[j] <- max(abs(check0.path.rml$est))
##   rvar0.path.rml[j] <- 
##   var(1/fips0.path.rml[iv==0,j])/mean(1/fips0.path.rml[iv==0,j])^2
## }
## par(mfrow=c(1,2))
## par(mar=c(4,4,2,2))
## plot(late.path.rml$ips[[2]]$nz.all, mdiff0.path.rml, 
##      xlim=c(0,p), ylim=c(0,.4), xlab="# nonzero", ylab="std diff")
## lines(late.path.rml$ips[[2]]$nz.all[!late.path.rml$ips[[2]]$non.conv],
##        mdiff0.path.rml, lty=3)
## points(late.path.rcal$ips[[1]]$nz.all[!late.path.rcal$ips[[1]]$non.conv],
##        mdiff0.path.rcal, pch=4)
## lines(late.path.rcal$ips[[1]]$nz.all[!late.path.rcal$ips[[1]]$non.conv],
##        mdiff0.path.rcal, lty=3)
## legend(120,.4, c("RML","RCAL"), pch=c(1,4), cex=.6)
## plot(rvar0.path.rml, mdiff0.path.rml, 
##      xlim=c(0,2.5), ylim=c(0,.4), xlab="rel var", ylab="std diff")
## lines(rvar0.path.rml, mdiff0.path.rml, lty=3)
## points(rvar0.path.rcal, mdiff0.path.rcal, pch=4)
## lines(rvar0.path.rcal, mdiff0.path.rcal, lty=3)
## legend(0.6,.4, c("RML","RCAL"), pch=c(1,4), cex=.6)


###################################################
### code chunk number 14: RCAL-LATE-vig.Rnw:491-514
###################################################
par(mfrow=c(1,2))
par(mar=c(4,4,2,2))
plot(late.path.rml$ips[[2]]$nz.all, mdiff1.path.rml, 
     xlim=c(0,p), ylim=c(0,.4), xlab="# nonzero", ylab="std diff")
lines(late.path.rml$ips[[2]]$nz.all[!late.path.rml$ips[[2]]$non.conv],
       mdiff1.path.rml, lty=3)

points(late.path.rcal$ips[[2]]$nz.all[!late.path.rcal$ips[[2]]$non.conv],
       mdiff1.path.rcal, pch=4)
lines(late.path.rcal$ips[[2]]$nz.all[!late.path.rcal$ips[[2]]$non.conv],
       mdiff1.path.rcal, lty=3)

legend(120,.4, c("RML","RCAL"), pch=c(1,4), cex=.6)

#
plot(rvar1.path.rml, mdiff1.path.rml, 
     xlim=c(0,1), ylim=c(0,.4), xlab="rel var", ylab="std diff")
lines(rvar1.path.rml, mdiff1.path.rml, lty=3)

points(rvar1.path.rcal, mdiff1.path.rcal, pch=4)
lines(rvar1.path.rcal, mdiff1.path.rcal, lty=3)

legend(0.6,.4, c("RML","RCAL"), pch=c(1,4), cex=.6)


