\name{simu.data}
\docType{data}
\alias{simu.data}

\title{Simulated data}
\description{A dataset simulated as in Tan (2019), Section 4.}

\usage{data(simu.data)}

\format{A data matrix with 800 rows and 202 columns.}

\details{ 
The dataset is generated as follows, where \code{y}, \code{tr}, and \code{x} represent an outcome, a treatment, and covariates respectively. 
\preformatted{
library(MASS)

###
mt0 <- 1-pnorm(-1)
mt1 <- dnorm(-1)
mt2 <- -(2*pnorm(-1)-1)/2 - dnorm(-1) +1/2
mt3 <- 3*dnorm(-1)
mt4 <- -3/2*(2*pnorm(-1)-1) - 4*dnorm(-1) +3/2

m.z1 <- mt0 + 2*mt1 + mt2
v.z1 <- mt0 + 4*mt1 + 6*mt2 + 4*mt3 + mt4
v.z1 <- v.z1 + 1 + 2*(mt1 + 2*mt2 + mt3)

sd.z1 <- sqrt(v.z1 -m.z1^2)
###

set.seed(123)

n <- 800
p <- 200

noise <- rnorm(n)

covm <- matrix(1,p,p)
for (i1 in 1:p)
  for (i2 in 1:p) {
    covm[i1,i2] <- 2^(-abs(i1-i2))
  }
x <- mvrnorm(n, mu=rep(0,p), Sigma=covm)

# transformation
z <- x 
for (i in 1:4) {
  z[,i] <- ifelse(x[,i]>-1,x[,i]+(x[,i]+1)^2,x[,i])
  z[,i] <- (z[,i]-m.z1) /sd.z1   # standardized
}

# treatment 
eta <- 1+ c( z[,1:4] \%*\% c(1, .5, .25, .125) )
tr <- rbinom(n, size=1, prob=expit(eta))

# outcome
eta.y <- c( z[,1:4] \%*\% c(1, .5, .25, .125) )
y <- eta.y + noise

# save 
simu.data <- cbind(y, tr, x)
save(simu.data, file="simu.data.rda")
} }

\references{

Tan, Z. (2019) Model-assisted inference for treatment effects using regularized calibrated estimation with high-dimensional data, \emph{Annals of Statistics}, to appear (preprint arXiv:1801.09817).
}

