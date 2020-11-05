\name{simu.iv.data}
\docType{data}
\alias{simu.iv.data}

\title{Simulated instrumental variable data}
\description{A dataset simulated as in Sun and Tan (2020), Section 4.}

\usage{data(simu.iv.data)}

\format{A data matrix with 800 rows and 203 columns.}

\details{ 
The dataset is generated as follows, where \code{y}, \code{iv}, \code{tr} and \code{x} represent an outcome, an instrumental variable, a treatment, and covariates respectively. 
\preformatted{
g<-function(z) {
	1/(1+exp(z/b))^2*dnorm(z)
}

rnorm.trunct <- function(n, mu, sig, lft, rgt) {
   x <- rep(0,n)
   for (i in 1:n) {
      x[i] <- rnorm(1,mu,sig)  
      while (x[i]<=lft | x[i]>rgt)
         x[i] <- rnorm(1,mu,sig)                                  
   }
   return(x)
}

### covariate mean and variance computed as in preprint of Tan (2020)

a<- 2.5;
c<- 2*pnorm(a)-1;
b<- sqrt(1-2*a*dnorm(a)/c)

m1<- exp(1/(8*b^2))*(pnorm(a-1/(2*b))-pnorm(-a-1/(2*b)))/c
v1<- exp(1/(2*b^2))*(pnorm(a-1/b)-pnorm(-a-1/b))/c-m1^2;

m2<- 10;
v2<- 1/c*integrate(g,-a,a)$value #by numerical integration

m3  <- 3/(25^2)*0.6+(0.6)^3;
mu4 <-(1/(b^4*c))*((3/2*(2*pnorm(a)-1)-a*(a^2+3)*dnorm(a))
	-(3/2*(2*pnorm(-a)-1)-(-a)*((-a)^2+3)*dnorm(-a)))
mu6 <-(1/(b^6*c))*((15/2*(2*pnorm(a)-1)-a*(a^4+5*a^2+15)*dnorm(a))
	-(15/2*(2*pnorm(-a)-1)-(-a)*((-a)^4+5*(-a)^2+15)*dnorm(-a)))
v3  <-mu6^2/25^6+15*mu4^2/25^4*0.6^2+15/25^2*0.6^4+0.6^6-m3^2

m4<- 2+20^2;
v4<- (2*mu4+6)+6*2*20^2+20^4-m4^2

###

set.seed(120)

n<- 800
p<- 200 

# covariates
 
x<- matrix(rnorm.trunct(p*n, 0, 1, -a, a),n,p)/b

# transformation

z<- x
z[,1] <- (exp(0.5*x[,1])-m1)/sqrt(v1);	
z[,2] <- (10+x[,2]/(1+exp(x[,1]))-m2)/sqrt(v2);	
z[,3] <- ((0.04*x[,1]*x[,3]+0.6)^3-m3)/sqrt(v3);
z[,4] <- ((x[,2]+x[,4]+20)^2-m4)/sqrt(v4);

# instrumental variable

eta<- z[,1:4]%*%c(1,-0.5,0.25,0.1)
iv<- rbinom(n,1,prob=expit(eta));

# unmeasured confounder in latent index model
u<- rlogis(n, location = 0, scale = 1);

# treatment
eta.d<- 1+cbind(iv,z[,1:4])%*%c(-2.5,0.25,1,0.5,-1.5)
tr<- as.numeric(eta.d >=u);

# outcome
late  <- 1
eta.y <- late*tr +z[,1:4]%*%c(0.5,1,1,1)+2*u
y <- rnorm(n, mean=eta.y, sd=1)

# save; if using main effects of x, then both the instrument propensity score
# and outcome models are misspecified

simu.iv.data <- cbind(y,tr,iv,x)

save(simu.iv.data, file="simu.iv.data.rda")
} }

\references{

Tan, Z. (2020) Model-assisted inference for treatment effects using regularized calibrated estimation with high-dimensional data, \emph{Annals of Statistics}, 48, 811–837.

Sun, B. and Tan, Z. (2020) High-dimensional model-assisted inference for local average treatment effects with instrumental variables, arXiv:2009.09286.

}

