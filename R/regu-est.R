
#' @import stats trust
NULL

### common functions
# expit()
# my.signif()
# n.signif()

# mn.ipw()
# mn.aipw()
# ate.ipw()
# ate.aipw()

### model-assisted inference with cross validation or along a regularization path
# mn.regu.cv()
# ate.regu.cv()

# mn.regu.path()
# ate.regu.path()

### model-assisted inference with non-regularized estimation
# mn.nreg()
# ate.nreg()

### non-regularized estimation using trust 
# loss.cal()
# loss.bal()
# glm.nreg()

### regularized estimation with cross validation or along a regularization path
# glm.regu.cv()
# glm.regu.path()

### linear and logistic (ml and cal) Lasso, also handling intercept and test data
# glm.regu()

### linear Lasso
# shift.down()
# shift.up()
# patch()
# asd()

###################################################################################
### common functions
###################################################################################

expit <- function(x)
  1/(1+exp(-x))

# ensure rounding up (not down)
#> signif(12.346, digits=3)
#[1] 12.3
#> my.signif(12.346, digits=3)
#[1] 12.4
#> signif(.12346, digits=3)
#[1] 0.123
#> my.signif(.12346, digits=3)
#[1] 0.124

my.signif <- function(x, digits=3) {

  y <- signif(x, digits)
  ifelse(x<=y, y, y +10^n.signif(abs(x-y)))
}

n.signif <- function(x) {
  # must x>0 

  if (x >=1) {
   i <- 1
   while (x>= 10^i) 
    i <- i+1
  } else {
   i <- 0
   while (x< 10^(i-1)) 
    i <- i-1
  }

  i
}

# IPW (ratio) estimator 

#' Inverse probability weighted estimation of population means
#'
#' This function implements inverse probability weighted (IPW) estimation of population means with missing data,
#' provided fitted propensity scores.
#'
#' The ratio IPW estimate is the direct IPW estimate divided by that with \code{y} replaced by a vector of 1s. The latter is referred to as
#' the direct IPW estimate of 1.  
#'
#' @param y An \eqn{n} x \eqn{1} vector of outcomes with missing data.
#' @param tr An \eqn{n} x \eqn{1} vector of non-missing indicators (=1 if \code{y} is observed or 0 if \code{y} is missing).
#' @param fp An \eqn{n} x \eqn{1} vector of fitted propensity scores.
#'
#' @return
#' \item{one}{The direct IPW estimate of 1.}
#' \item{est}{The ratio IPW estimate.}
#'
#' @references
#' Tan, Z. (2017) Regularized calibrated estimation of propensity scores with model misspecification and high-dimensional data, arXiv:1710.08074. 
#'
#' Tan, Z. (2019) Model-assisted inference for treatment effects using regularized calibrated estimation with high-dimensional data, 
#' \emph{Annals of Statistics}, to appear (preprint arXiv:1801.09817).
#'
#' @export

mn.ipw <- function(y, tr, fp) {

  y <- cbind(y)

  n <- length(tr)
  one <- sum(1/fp[tr==1]) /n

  est <- apply(y[tr==1, ,drop=F] /fp[tr==1], 2, sum) /n

  list(one=one, est=est /one)
}


# Augmented IPW estimator 

#' Augmented inverse probability weighted estimation of population means
#'
#' This function implements augmented inverse probability weighted (IPW) estimation of population means with missing data,
#' provided both fitted propensity scores and fitted values from outcome regression.
#'
#' @param y An \eqn{n} x \eqn{1} vector of outcomes with missing data.
#' @param tr An \eqn{n} x \eqn{1} vector of non-missing indicators (=1 if \code{y} is observed or 0 if \code{y} is missing).
#' @param fp An \eqn{n} x \eqn{1} vector of fitted propensity scores.
#' @param fo An \eqn{n} x \eqn{1} vector of fitted values from outcome regression.
#' @param off An offset value (e.g., the true value in simulations) used to calculate the z-statistic.
#'
#' @return
#' \item{one}{The direct IPW estimate of 1.}
#' \item{ipw}{The ratio IPW estimate.}
#' \item{or}{The outcome regression estimate.}
#' \item{est}{The augmented IPW estimate.}
#' \item{var}{The estimated variance associated with the augmented IPW estimate.}
#' \item{ze}{The z-statistic for the augmented IPW estimate, compared to \code{off}.}
#'
#' @references
#' Tan, Z. (2017) Regularized calibrated estimation of propensity scores with model misspecification and high-dimensional data, arXiv:1710.08074. 
#'
#' Tan, Z. (2019) Model-assisted inference for treatment effects using regularized calibrated estimation with high-dimensional data, 
#' \emph{Annals of Statistics}, to appear (preprint arXiv:1801.09817).
#'
#' @export

mn.aipw <- function(y, tr, fp, fo, off=0) {

  n <- length(tr)
  one <- sum(1/fp[tr==1]) /n

  est1 <- sum(y[tr==1] /fp[tr==1]) /n

  est2 <- mean(fo)

  est3 <- sum((y[tr==1]-fo[tr==1]) /fp[tr==1]) /n + mean(fo)
  var3 <- mean( ifelse(tr==1, (y-fo)/fp +fo-est3, fo-est3)^2 ) /n

  list(one=one, ipw=est1 /one, or=est2, 
       est=est3, var=var3, ze=(est3-off)/sqrt(var3))
}

#' Inverse probability weighted estimation of average treatment effects
#'
#' This function implements inverse probability weighted (IPW) estimation of average treatment effects (ATEs),
#' provided fitted propensity scores.
#'
#' @param y An \eqn{n} x \eqn{1} vector of observed outcomes.
#' @param tr An \eqn{n} x \eqn{1} vector of treatment indicators (=1 if treated or 0 if untreated).
#' @param mfp An \eqn{n} x \eqn{2} matrix of fitted propensity scores for untreated (first column) and treated (second column).
#'
#' @return
#' \item{one}{The direct IPW estimates of 1.}
#' \item{est}{The ratio IPW estimates of means.}
#' \item{diff}{The ratio IPW estimate of ATE.}
#'
#' @references
#' Tan, Z. (2017) Regularized calibrated estimation of propensity scores with model misspecification and high-dimensional data, arXiv:1710.08074. 
#'
#' Tan, Z. (2019) Model-assisted inference for treatment effects using regularized calibrated estimation with high-dimensional data, 
#' \emph{Annals of Statistics}, to appear (preprint arXiv:1801.09817).
#'
#' @export

ate.ipw <- function(y, tr, mfp) {

  #d takes value 1,2,... 
  d <- 1+tr

  n <- length(d)
  m <- dim(mfp)[2]

  one <- rep(NA, m)
  est <- rep(NA, m)
  diff.est <- rep(NA, m)

  for (j in 1:dim(mfp)[2]) {

    one[j] <- sum(1/mfp[d==j,j]) /n

    est[j] <- sum(y[d==j] /mfp[d==j,j]) /n /one[j]

    if (j>1) 
      diff.est[j] <- est[j] -est[1]
  }

  list(one=one, est=est, diff.est=diff.est)
}


#' Augmented inverse probability weighted estimation of population means
#'
#' This function implements augmented inverse probability weighted (IPW) estimation of average treatment effects (ATEs),
#' provided both fitted propensity scores and fitted values from outcome regression.
#'
#' @param y An \eqn{n} x \eqn{1} vector of observed outcomes.
#' @param tr An \eqn{n} x \eqn{1} vector of treatment indicators (=1 if treated or 0 if untreated).
#' @param mfp An \eqn{n} x \eqn{2} matrix of fitted propensity scores for untreated (first column) and treated (second column).
#' @param mfo An \eqn{n} x \eqn{2} matrix of fitted values from outcome regression, for untreated (first column) and treated (second column).
#' @param off A \eqn{2} x \eqn{1} vector of offset values (e.g., the true values in simulations) used to calculate the z-statistics.
#'
#' @return
#' \item{one}{A \eqn{2} x \eqn{1} vector of direct IPW estimates of 1.}
#' \item{ipw}{A \eqn{2} x \eqn{1} vector of ratio IPW estimates of means.}
#' \item{or}{A \eqn{2} x \eqn{1} vector of outcome regression estimates of means.}
#' \item{est}{A \eqn{2} x \eqn{1} vector of augmented IPW estimates of means.}
#' \item{var}{The estimated variances associated with the augmented IPW estimates of means.}
#' \item{ze}{The z-statistics for the augmented IPW estimates of means, compared to \code{off}.}
#' \item{diff}{The augmented IPW estimate of ATE.}
#' \item{diff.var}{The estimated variance associated with the augmented IPW estimate of ATE.}
#' \item{diff.ze}{The z-statistic for the augmented IPW estimate of ATE.}
#'
#' @references
#' Tan, Z. (2017) Regularized calibrated estimation of propensity scores with model misspecification and high-dimensional data, arXiv:1710.08074. 
#'
#' Tan, Z. (2019) Model-assisted inference for treatment effects using regularized calibrated estimation with high-dimensional data, 
#' \emph{Annals of Statistics}, to appear (preprint arXiv:1801.09817).
#'
#' @export

ate.aipw <- function(y, tr, mfp, mfo, off=NULL) {

  #d takes value 1,2,... 
  d <- 1+tr

  n <- length(d)
  m <- dim(mfp)[2]

  if (is.null(off)) off <- rep(0,m)

  one <- rep(NA, m)
  est1 <- rep(NA, m)
  est2 <- rep(NA, m)
  est3 <- rep(NA, m)
  var3 <- rep(NA, m)

  phi <- matrix(0, nrow=n, ncol=m)

  est4 <- rep(NA, m)
  var4 <- rep(NA, m)  

  for (j in 1:dim(mfp)[2]) {

    one[j] <- sum(1/mfp[d==j,j]) /n

    est1[j] <- sum(y[d==j] /mfp[d==j,j]) /n

    est2[j] <- mean(mfo[,j])

    phi[,j] <- mfo[,j]
    phi[d==j,j] <- phi[d==j,j] + (y[d==j]-mfo[d==j,j]) /mfp[d==j,j]

    est3[j] <- mean(phi[,j])

    var3[j] <- mean( (phi[,j]-est3[j])^2 ) /n

    if (j>1) {
      est4[j] <- est3[j] - est3[1]
      var4[j] <- mean( (phi[,j]-phi[,1] -(est3[j]-est3[1]))^2 ) /n
    }
  }

  list(one=one, ipw=est1 /one, or=est2, 
       est=est3, var=var3, ze=(est3-off)/sqrt(var3),
       diff.est=est4, diff.var=var4, diff.ze=(est4-(off-off[1]))/sqrt(var4))
}

###################################################################################
### model-assisted inference with cross validation
###################################################################################

#' Model-assisted inference for population means based on cross validation
#'
#' This function implements model-assisted inference for population means with missing data,
#' using regularized calibrated estimation based on cross validation.
#'
#' Two steps are involved in this function: first fitting propensity score and outcome regression models and then applying the augmented IPW estimator 
#' for a population mean. For \code{ploss}="cal", regularized calibrated estimation is performed with cross validation as described in Tan (2017, 2019). 
#' The method then leads to model-assisted inference, in which confidence intervals are valid with high-dimensinoal data 
#' if the propensity score model is correctly specified but the outcome regression model may be misspecified.
#' With linear outcome models, the inference is also doubly robust.
#' For \code{ploss}="ml", regularized maximum likelihood estimation is used (Belloni et al. 2014; Farrell 2015). In this case, standard errors 
#' are only shown to be valid if both the propensity score model and the outcome regression model are correctly specified.
#'
#' @param fold A vector of length 2 giving the fold numbers for cross validation in propensity score estimation and outcome regression respectively.
#' @param nrho A vector of length 2 giving the numbers of tuning parameters searched in cross validation.
#' @param rho.seq A list of two vectors giving the tuning parameters in propensity score estimation (first vector) and outcome regression (second vector). 
#' @param y An \eqn{n} x \eqn{1} vector of outcomes with missing data.
#' @param tr An \eqn{n} x \eqn{1} vector of non-missing indicators (=1 if \code{y} is observed or 0 if \code{y} is missing).
#' @param x An \eqn{n} x \eqn{p} matix of covariates, used in both propensity score and outcome regression models. 
#' @param ploss A loss function used in propensity score estimation (either "ml" or "cal").
#' @param yloss A loss function used in outcome regression (either "gaus" for continuous outcomes or "ml" for binary outcomes).
#' @param off An offset value (e.g., the true value in simulations) used to calculate the z-statistic from augmented IPW estimation.
#' @param ... Additional arguments to \code{\link{glm.regu.cv}}.
#'
#' @return
#' \item{ps}{A list containing the results from fitting the propensity score model by \code{\link{glm.regu.cv}}.}
#' \item{fp}{The \eqn{n} x \eqn{1} vector of fitted propensity scores.}
#' \item{or}{A list containing the results from fitting the outcome regression model by \code{\link{glm.regu.cv}}.}
#' \item{fo}{The \eqn{n} x \eqn{1} vector of fitted values from outcome regression.}
#' \item{est}{A list containing the results from augmented IPW estimation by \code{\link{mn.aipw}}.}
#'
#' @references
#' Belloni, A., Chernozhukov, V., and Hansen, C. (2014) Inference on treatment effects after selection among high-dimensional controls,
#' \emph{Review of Economic Studies}, 81, 608-650.
#'
#' Farrell, M.H. (2015) Robust inference on average treatment effects with possibly more covariates than observations, \emph{Journal of Econometrics}, 189, 1-23.
#'
#' Tan, Z. (2017) Regularized calibrated estimation of propensity scores with model misspecification and high-dimensional data, arXiv:1710.08074. 
#'
#' Tan, Z. (2019) Model-assisted inference for treatment effects using regularized calibrated estimation with high-dimensional data, 
#' \emph{Annals of Statistics}, to appear (preprint arXiv:1801.09817).
#'
#' @examples 
#' \donttest{
#' data(simu.data)
#' n <- dim(simu.data)[1]
#' p <- dim(simu.data)[2]-2
#'
#' y <- simu.data[,1]
#' tr <- simu.data[,2]
#' x <- simu.data[,2+1:p]
#' x <- scale(x)
#'
#' # missing data
#' y[tr==0] <- NA
#'
#' mn.cv.rcal <- mn.regu.cv(fold=5*c(1,1), nrho=(1+10)*c(1,1), rho.seq=NULL, y, tr, x, 
#'                          ploss="cal", yloss="gaus")
#' unlist(mn.cv.rcal$est)
#' }
#'
#' @export

mn.regu.cv <- function(fold, nrho=NULL, rho.seq=NULL, y, tr, x,
                       ploss="cal", yloss="gaus", off=0, ...) {

  # propensity score
  ps.regu <- glm.regu.cv(fold=fold[1], nrho=nrho[1], rho.seq=rho.seq[[1]], y=tr, x=x, 
                         iw=NULL, loss=ploss, ...)

  fp <- ps.regu $sel.fit[,1]
  
  # ourcome regression
  if (ploss=="ml") {
    iw <- NULL
  } else {  # "cal"
    iw <- 1/fp[tr==1] -1
  }

  or.regu <- glm.regu.cv(fold=fold[2], nrho=nrho[2], rho.seq=rho.seq[2], y[tr==1], x[tr==1,], 
                         iw=iw, loss=yloss, ...)

  if (yloss=="gaus") {
    fo <- c( cbind(1,x)%*%or.regu$sel.bet[,1] )
  } else { # binary y
    fo <- expit( c( cbind(1,x)%*%or.regu$sel.bet[,1] ) )
  }

  # AIPW estimation
  out.aipw <- mn.aipw(y, tr, fp=fp, fo=fo, off=off)

  list(ps=ps.regu, fp=fp, or=or.regu, fo=fo, est=out.aipw)
}

#' Model-assisted inference for average treatment effects based on cross validation
#'
#' This function implements model-assisted inference for average treatment effects,
#' using regularized calibrated estimation based on cross validation.
#'
#' For calibrated estimation, two sets of propensity scores are separately estimated for the untreated and treated as discussed in Tan (2017, 2019).
#' See also \strong{Details} for \code{\link{mn.regu.cv}}.
#'
#' @param fold A vector of length 2 giving the fold numbers for cross validation in propensity score estimation and outcome regression respectively.
#' @param nrho A vector of length 2 giving the numbers of tuning parameters searched in cross validation.
#' @param rho.seq A list of two vectors giving the tuning parameters in propensity score estimation (first vector) and outcome regression (second vector). 
#' @param y An \eqn{n} x \eqn{1} vector of obseved outcomes.
#' @param tr An \eqn{n} x \eqn{1} vector of treatment indicators (=1 if treated or 0 if untreated).
#' @param x An \eqn{n} x \eqn{p} matix of covariates, used in both propensity score and outcome regression models. 
#' @param ploss A loss function used in propensity score estimation (either "ml" or "cal").
#' @param yloss A loss function used in outcome regression (either "gaus" for continuous outcomes or "ml" for binary outcomes).
#' @param off A \eqn{2} x \eqn{1} vector of offset values (e.g., the true values in simulations) used to calculate the z-statistics from augmented IPW estimation.
#' @param ... Additional arguments to \code{\link{glm.regu.cv}}.
#'
#' @return
#' \item{ps}{A list containing the results from fitting the propensity score model by \code{\link{glm.regu.cv}}.}
#' \item{mfp}{An \eqn{n} x \eqn{2} matrix of fitted propensity scores for untreated (first column) and treated (second column).}
#' \item{or}{A list containing the results from fitting the outcome regression model by \code{\link{glm.regu.cv}}.}
#' \item{mfo}{An \eqn{n} x \eqn{2} matrix of fitted values from outcome regression, for untreated (first column) and treated (second column).}
#' \item{est}{A list containing the results from augmented IPW estimation by \code{\link{ate.aipw}}.}
#'
#' @references
#' Tan, Z. (2017) Regularized calibrated estimation of propensity scores with model misspecification and high-dimensional data, arXiv:1710.08074. 
#'
#' Tan, Z. (2019) Model-assisted inference for treatment effects using regularized calibrated estimation with high-dimensional data, 
#' \emph{Annals of Statistics}, to appear (preprint arXiv:1801.09817).
#'
#' @examples 
#' \donttest{
#' data(simu.data)
#' n <- dim(simu.data)[1]
#' p <- dim(simu.data)[2]-2
#'
#' y <- simu.data[,1]
#' tr <- simu.data[,2]
#' x <- simu.data[,2+1:p]
#' x <- scale(x)
#'
#' ate.cv.rcal <- ate.regu.cv(fold=5*c(1,1), nrho=(1+10)*c(1,1), rho.seq=NULL, y, tr, x, 
#'                            ploss="cal", yloss="gaus")
#'
#' matrix(unlist(ate.cv.rcal$est), ncol=2, byrow=TRUE, 
#' dimnames=list(c("one", "ipw", "or", "est", "var", "ze", 
#' "diff.est", "diff.var", "diff.ze"), c("untreated", "treated")))
#' }
#'
#' @export

ate.regu.cv <- function(fold, nrho=NULL, rho.seq=NULL, y, tr, x,
                        ploss="cal", yloss="gaus", off=NULL, ...) {

  n <- length(tr)
  d <- 1+tr

  ps.regu <- vector("list", 2)
  or.regu <- vector("list", 2)

  mfp <- matrix(NA, nrow=n, ncol=2)
  mfo <- matrix(NA, nrow=n, ncol=2)

  for (k in 2:1) {

    # propensity score
    if (k==2 | ploss=="cal") {
      ps.regu[[k]] <- glm.regu.cv(fold=fold[1], nrho=nrho[1], rho.seq=rho.seq[[1]], y= d==k, x=x, 
                                  iw=NULL, loss=ploss, ...)

      mfp[,k] <- ps.regu[[k]] $sel.fit[,1]

    } else { # k==1 & ploss=="ml"

      mfp[,1] <- 1-mfp[,2]
    }

    # outcome regression
    if (ploss=="ml") {
      iw <- NULL
    } else {  # "cal"
      iw <- 1/mfp[d==k,k] -1
    }

    or.regu[[k]] <- glm.regu.cv(fold=fold[2], nrho=nrho[2], rho.seq=rho.seq[2], y[d==k], x[d==k,], 
                                iw=iw, loss=yloss, ...)

    if (yloss=="gaus") {
      mfo[,k] <- c( cbind(1,x)%*%or.regu[[k]] $sel.bet[,1] )
    } else { # binary y
      mfo[,k] <- expit( c( cbind(1,x)%*%or.regu[[k]] $sel.bet[,1] ) )
    }
  }  

  # AIPW estimation
  out.aipw <- ate.aipw(y, tr, mfp=mfp, mfo=mfo, off=off)

  list(ps=ps.regu, mfp=mfp, or=or.regu, mfo=mfo, est=out.aipw)
}


###################################################################################
### model-assisted inference along a regularization path for PS (cross validation for OR)
###################################################################################

#' Model-assisted inference for population means along a regularization path
#'
#' This function implements model-assisted inference for population means with missing data,
#' using regularized calibrated estimation along a regularization path for propensity score (PS) estimation 
#' while based on cross validation for outcome regression (OR).
#'
#' See \strong{Details} for \code{\link{mn.regu.cv}}.
#'
#' @param fold A vector of length 2, with the second component giving the fold number for cross validation in outcome regression. The first component is not used.
#' @param nrho A vector of length 2 giving the number of tuning parameters in a regularization path for PS estimation and that in cross validation for OR.
#' @param rho.seq A list of two vectors giving the tuning parameters for propensity score estimation (first vector) and outcome regression (second vector). 
#' @param y An \eqn{n} x \eqn{1} vector of outcomes with missing data.
#' @param tr An \eqn{n} x \eqn{1} vector of non-missing indicators (=1 if \code{y} is observed or 0 if \code{y} is missing).
#' @param x An \eqn{n} x \eqn{p} matix of covariates, used in both propensity score and outcome regression models. 
#' @param ploss A loss function used in propensity score estimation (either "ml" or "cal").
#' @param yloss A loss function used in outcome regression (either "gaus" for continuous outcomes or "ml" for binary outcomes).
#' @param off An offset value (e.g., the true value in simulations) used to calculate the z-statistic from augmented IPW estimation.
#' @param ... Additional arguments to \code{\link{glm.regu.cv}} and \code{\link{glm.regu.path}}.
#'
#' @return
#' \item{ps}{A list containing the results from fitting the propensity score model by \code{\link{glm.regu.path}}.}
#' \item{fp}{The matrix of fitted propensity scores, column by column, along the PS regularization path.}
#' \item{or}{A list of objects, each giving the results from fitting the outcome regression model by \code{\link{glm.regu.cv}} for a PS tuning parameter.}
#' \item{fo}{The matrix of fitted values from outcome regression based on cross validation, column by column, along the PS regularization path.}
#' \item{est}{A list containing the results from augmented IPW estimation by \code{\link{mn.aipw}}.}
#' \item{rho}{A vector of tuning parameters leading to converged results in propensity score estimation.}
#'
#' @references
#' Tan, Z. (2017) Regularized calibrated estimation of propensity scores with model misspecification and high-dimensional data, arXiv:1710.08074. 
#'
#' Tan, Z. (2019) Model-assisted inference for treatment effects using regularized calibrated estimation with high-dimensional data, 
#' \emph{Annals of Statistics}, to appear (preprint arXiv:1801.09817).
#'
#' @examples 
#' \donttest{
#' data(simu.data)
#' n <- dim(simu.data)[1]
#' p <- dim(simu.data)[2]-2
#'
#' y <- simu.data[,1]
#' tr <- simu.data[,2]
#' x <- simu.data[,2+1:p]
#' x <- scale(x)
#'
#' # missing data
#' y[tr==0] <- NA
#'
#' mn.path.rcal <- mn.regu.path(fold=5*c(0,1), nrho=(1+10)*c(1,1), y=y, tr=tr, x=x, 
#'                              ploss="cal", yloss="gaus")
#' mn.path.rcal$est
#' }
#' 
#' @export

mn.regu.path <- function(fold, nrho=NULL, rho.seq=NULL, y, tr, x,
                         ploss="cal", yloss="gaus", off=0, ...) {

  n <- length(tr)

  # propensity score
  ps.regu <- glm.regu.path(nrho=nrho[1], rho.seq=rho.seq[[1]], y=tr, x=x, 
                           iw=NULL, loss=ploss, ...)

  fp.all <- ps.regu $fit.all[, !ps.regu$non.conv, drop=F]

  # outcome regression
  or.regu <- vector("list", dim(fp.all)[2])
  fo.all <- matrix(NA, nrow=n, ncol=dim(fp.all)[2])

  est.all <- matrix(NA, nrow=6, ncol=dim(fp.all)[2])
  rownames(est.all) <- c("one", "ipw", "or", "est", "var", "ze")

  for (j in 1:dim(fp.all)[2]) {

    if (ploss=="ml") {
      iw <- NULL
      if (j==1) {
        or.regu[[j]] <- glm.regu.cv(fold=fold[2], nrho=nrho[2], rho.seq=rho.seq[2], y[tr==1], x[tr==1,], 
                                    iw=iw, loss=yloss, ...)
      } else {
        or.regu[[j]] <- or.regu[[1]]
      }

    } else {  # "cal"
      iw <- 1/fp.all[tr==1,j] -1
      or.regu[[j]] <- glm.regu.cv(fold=fold[2], nrho=nrho[2], rho.seq=rho.seq[2], y[tr==1], x[tr==1,], 
                                  iw=iw, loss=yloss, ...)
    }


    if (yloss=="gaus") {
      fo.all[,j] <- c( cbind(1,x)%*%or.regu[[j]]$sel.bet[,1] )
    } else { # binary y
      fo.all[,j] <- expit( c( cbind(1,x)%*%or.regu[[j]]$sel.bet[,1] ) )
    }

    # AIPW estimation
    est.all[,j] <- unlist(mn.aipw(y, tr, fp=fp.all[,j], fo=fo.all[,j], off=off))
  }

  list(ps=ps.regu, fp=fp.all, or=or.regu, fo=fo.all, est=est.all, rho=ps.regu$rho[!ps.regu$non.conv])
}

#' Model-assisted inference for average treatment effects along regularization paths
#'
#' This function implements model-assisted inference for average treatment effects,
#' using regularized calibrated estimation along regularization paths for propensity score (PS) estimation 
#' while based on cross validation for outcome regression (OR).
#'
#' See \strong{Details} for \code{\link{ate.regu.cv}}.
#'
#' @param fold A vector of length 2, with the second component giving the fold number for cross validation in outcome regression. The first component is not used.
#' @param nrho A vector of length 2 giving the number of tuning parameters in a regularization path for PS estimation and that in cross validation for OR.
#' @param rho.seq A list of two vectors giving the tuning parameters for propensity score estimation (first vector) and outcome regression (second vector). 
#' @param y An \eqn{n} x \eqn{1} vector of observed outcomes.
#' @param tr An \eqn{n} x \eqn{1} vector of treatment indicators (=1 if treated or 0 if untreated).
#' @param x An \eqn{n} x \eqn{p} matix of covariates, used in both propensity score and outcome regression models. 
#' @param ploss A loss function used in propensity score estimation (either "ml" or "cal").
#' @param yloss A loss function used in outcome regression (either "gaus" for continuous outcomes or "ml" for binary outcomes).
#' @param off A \eqn{2} x \eqn{1} vector of offset values (e.g., the true values in simulations) used to calculate the z-statistics from augmented IPW estimation.
#' @param ... Additional arguments to \code{\link{glm.regu.cv}} and \code{\link{glm.regu.path}}.
#'
#' @return
#' \item{ps}{A list of 2 objects, giving the results from fitting the propensity score model by \code{\link{glm.regu.path}} for untreated (first) and treated (second).}
#' \item{mfp}{A list of 2 matrices of fitted propensity scores, along the PS regularization path, for untreated (first matrix) and treated (second matrix).}
#' \item{or}{A list of 2 lists of objects for untreated (first) and treated (second), where each object gives 
#' the results from fitting the outcome regression model by \code{\link{glm.regu.cv}} for a PS tuning parameter.}
#' \item{mfo}{A list of 2 matrices of fitted values from outcome regression based on cross validation, along the PS regularization path, 
#' for untreated (first matrix) and treated (second matrix).}
#' \item{est}{A list containing the results from augmented IPW estimation by \code{\link{ate.aipw}}.}
#' \item{rho}{A vector of tuning parameters leading to converged results in propensity score estimation.}
#'
#' @references
#' Tan, Z. (2017) Regularized calibrated estimation of propensity scores with model misspecification and high-dimensional data, arXiv:1710.08074. 
#'
#' Tan, Z. (2019) Model-assisted inference for treatment effects using regularized calibrated estimation with high-dimensional data, 
#' \emph{Annals of Statistics}, to appear (preprint arXiv:1801.09817).
#'
#' @examples 
#' \donttest{
#' data(simu.data)
#' n <- dim(simu.data)[1]
#' p <- dim(simu.data)[2]-2
#'
#' y <- simu.data[,1]
#' tr <- simu.data[,2]
#' x <- simu.data[,2+1:p]
#' x <- scale(x)
#'
#' ate.path.rcal <- ate.regu.path(fold=5*c(0,1), nrho=(1+10)*c(1,1), rho.seq=NULL, y, tr, x, 
#'                                ploss="cal", yloss="gaus")
#' ate.path.rcal$est
#' }
#'
#' @export

ate.regu.path <- function(fold, nrho=NULL, rho.seq=NULL, y, tr, x,
                          ploss="cal", yloss="gaus", off=NULL, ...) {

  n <- length(tr)
  d <- 1+tr

  ps.regu <- vector("list", 2)
  or.regu <- vector("list", 2)

  mfp.all <- vector("list", 2)
  mfo.all <- vector("list", 2)

  for (k in 2:1) {

    # propensity score
    if (k==2 | ploss=="cal") {
      ps.regu[[k]] <- glm.regu.path(nrho=nrho[1], rho.seq=rho.seq[[1]], y= d==k, x=x, 
                                    iw=NULL, loss=ploss, ...)

      mfp.all[[k]] <- ps.regu[[k]] $fit.all[, !ps.regu[[k]]$non.conv, drop=F]

    } else { # k==1 & ploss=="ml"

      mfp.all[[1]] <- 1-mfp.all[[2]]
    }

    # ourcome regression
    or.regu[[k]] <- vector("list", dim(mfp.all[[k]])[2])
    mfo.all[[k]] <- matrix(NA, nrow=n, ncol=dim(mfp.all[[k]])[2])

    for (j in 1:dim(mfp.all[[k]])[2]) {
      if (ploss=="ml") {
        iw <- NULL
        if (j==1) {
          or.regu[[k]][[j]] <- glm.regu.cv(fold=fold[2], nrho=nrho[2], rho.seq=rho.seq[2], y[d==k], x[d==k,], 
                                           iw=iw, loss=yloss, ...)
        } else {
          or.regu[[k]][[j]] <- or.regu[[k]][[1]]
        }

      } else {  # "cal"
        iw <- 1/mfp.all[[k]][d==k,j] -1
        or.regu[[k]][[j]] <- glm.regu.cv(fold=fold[2], nrho=nrho[2], rho.seq=rho.seq[2], y[d==k], x[d==k,], 
                                         iw=iw, loss=yloss, ...)
      }

      if (yloss=="gaus") {
        mfo.all[[k]][,j] <- c( cbind(1,x)%*%or.regu[[k]][[j]] $sel.bet[,1] )
      } else { # binary y
        mfo.all[[k]][,j] <- expit( c( cbind(1,x)%*%or.regu[[k]][[j]] $sel.bet.all[,1] ) )
      }
    }
  }  

  # AIPW estimation
  num.conv <- c(dim(mfp.all[[1]])[2], dim(mfp.all[[2]])[2])
  est.all <- array(NA, c(9, 2, min(num.conv)))
  dimnames(est.all)[[1]] <- c("one", "ipw", "or", "est", "var", "ze", "diff.est", "diff.var", "diff.ze")

  for (j in 1:min(num.conv)) {
     est.all[,, min(num.conv)+1-j] <- 
       matrix(unlist(ate.aipw(y, tr, mfp=cbind(mfp.all[[1]][,num.conv[1]+1-j],mfp.all[[2]][,num.conv[2]+1-j]), 
          mfo=cbind(mfo.all[[1]][,num.conv[1]+1-j],mfo.all[[2]][,num.conv[2]+1-j]), off=off)), ncol=2,byrow=TRUE)
  }

  list(ps=ps.regu, mfp=mfp.all, or=or.regu, mfo=mfo.all, est=est.all, 
       rho=ps.regu[[2]]$rho[nrho[1]-min(num.conv) +1:min(num.conv)]) #ps.regu$rho[[1]] is NULL if ploss=="ml"
}


###################################################################################
### model-assisted inference with non-regularized estimation
###################################################################################

#' Model-assisted inference for population means without regularization 
#'
#' This function implements model-assisted inference for population means with missing data,
#' using non-regularized calibrated estimation.
#'
#' Two steps are involved in this function: first fitting propensity score and outcome regression models and then applying the augmented IPW estimator 
#' for a population mean. For \code{ploss}="cal", calibrated estimation is performed similarly as in Tan (2017, 2019), but without regularization. 
#' The method then leads to model-assisted inference, in which confidence intervals are valid if the propensity score model is correctly specified but 
#' the outcome regression model may be misspecified.
#' With linear outcome models, the inference is also doubly robust (Kim and Haziza 2014; Vermeulen and Vansteelandt 2015).  
#' For \code{ploss}="ml", maximum likelihood estimation is used (Robins et al. 1994). In this case, standard errors are in general conservative 
#' if the propensity score model is correctly specified but the outcome regression model may be misspecified.
#'
#' @param y An \eqn{n} x \eqn{1} vector of outcomes with missing data.
#' @param tr An \eqn{n} x \eqn{1} vector of non-missing indicators (=1 if \code{y} is observed or 0 if \code{y} is missing).
#' @param x An \eqn{n} x \eqn{p} matix of covariates (excluding a constant), used in both propensity score and outcome regression models. 
#' @param ploss A loss function used in propensity score estimation (either "ml" or "cal").
#' @param yloss A loss function used in outcome regression (either "gaus" for continuous outcomes or "ml" for binary outcomes).
#' @param off An offset value (e.g., the true value in simulations) used to calculate the z-statistic from augmented IPW estimation.
#'
#' @return
#' \item{ps}{A list containing the results from fitting the propensity score model by \code{\link{glm.nreg}}.}
#' \item{fp}{The \eqn{n} x \eqn{1} vector of fitted propensity scores.}
#' \item{or}{A list containing the results from fitting the outcome regression model by \code{\link{glm.nreg}}.}
#' \item{fo}{The \eqn{n} x \eqn{1} vector of fitted values from outcome regression.}
#' \item{est}{A list containing the results from augmented IPW estimation by \code{\link{mn.aipw}}.}
#'
#' @references
#' Kim, J.K. and Haziza, D. (2014) Doubly robust inference with missing data in survey sampling, \emph{Statistica Sinica}, 24, 375-394.
#'
#' Robins, J.M., Rotnitzky, A., and Zhao, L.P. (1994) Estimation of regression coefficients when some regressors are not always observed, 
#' \emph{Journal of the American Statistical Association}, 89, 846-866.
#'
#' Tan, Z. (2017) Regularized calibrated estimation of propensity scores with model misspecification and high-dimensional data, arXiv:1710.08074. 
#'
#' Tan, Z. (2019) Model-assisted inference for treatment effects using regularized calibrated estimation with high-dimensional data, 
#' \emph{Annals of Statistics}, to appear (preprint arXiv:1801.09817).
#'
#' Vermeulen, K. and Vansteelandt, S. (2015) Bias-reduced doubly robust estimation, \emph{Journal of the American Statistical Association}, 110, 1024-1036.
#'
#' @examples 
#' data(simu.data)
#' n <- dim(simu.data)[1]
#' p <- dim(simu.data)[2]-2
#'
#' y <- simu.data[,1]
#' tr <- simu.data[,2]
#' x <- simu.data[,2+1:p]
#' x <- scale(x)
#'
#' # missing data
#' y[tr==0] <- NA
#'
#' # include only 10 covariates
#' x2 <- x[,1:10]
#'
#' mn.cal <- mn.nreg(y, tr, x2, ploss="cal", yloss="gaus")
#' unlist(mn.cal$est)
#'
#' @export

mn.nreg <- function(y, tr, x, ploss="cal", yloss="gaus", off=0) {

  # propensity score
  ps.nreg <- glm.nreg(y=tr, x=x, iw=NULL, loss=ploss)

  fp <- ps.nreg $fit
  
  # ourcome regression
  if (ploss=="ml") {
    iw <- NULL
  } else {  # "cal"
    iw <- 1/fp[tr==1] -1
  }

  or.nreg <- glm.nreg(y[tr==1], x[tr==1,], iw=iw, loss=yloss)

  if (yloss=="gaus") {
    fo <- c( cbind(1,x)%*%or.nreg$coef )
  } else { # binary y
    fo <- expit( c( cbind(1,x)%*%or.nreg$coef ) )
  }

  # AIPW estimation
  out.aipw <- mn.aipw(y, tr, fp=fp, fo=fo, off=off)

  list(ps=ps.nreg, fp=fp, or=or.nreg, fo=fo, est=out.aipw)
}

#' Model-assisted inference for average treatment effects without regularization 
#'
#' This function implements model-assisted inference for average treatment effects,
#' using non-regularized calibrated estimation.
#'
#' For calibrated estimation, two sets of propensity scores are separately estimated for the untreated and treated as discussed in Tan (2017, 2019).
#' See also \strong{Details} for \code{\link{mn.nreg}}.
#'
#' @param y An \eqn{n} x \eqn{1} vector of observed outcomes.
#' @param tr An \eqn{n} x \eqn{1} vector of treatment indicators (=1 if treated or 0 if untreated).
#' @param x An \eqn{n} x \eqn{p} matix of covariates, used in both propensity score and outcome regression models. 
#' @param ploss A loss function used in propensity score estimation (either "ml" or "cal").
#' @param yloss A loss function used in outcome regression (either "gaus" for continuous outcomes or "ml" for binary outcomes).
#' @param off A \eqn{2} x \eqn{1} vector of offset values (e.g., the true values in simulations) used to calculate the z-statistics from augmented IPW estimation.
#'
#' @return
#' \item{ps}{A list containing the results from fitting the propensity score model by \code{\link{glm.nreg}}.}
#' \item{mfp}{An \eqn{n} x \eqn{2} matrix of fitted propensity scores for untreated (first column) and treated (second column).}
#' \item{or}{A list containing the results from fitting the outcome regression model by \code{\link{glm.nreg}}.}
#' \item{mfo}{An \eqn{n} x \eqn{2} matrix of fitted values from outcome regression, for untreated (first column) and treated (second column).}
#' \item{est}{A list containing the results from augmented IPW estimation by \code{\link{ate.aipw}}.}
#'
#' @references
#' Tan, Z. (2017) Regularized calibrated estimation of propensity scores with model misspecification and high-dimensional data, arXiv:1710.08074. 
#'
#' Tan, Z. (2019) Model-assisted inference for treatment effects using regularized calibrated estimation with high-dimensional data, 
#' \emph{Annals of Statistics}, to appear (preprint arXiv:1801.09817).
#'
#' @examples 
#' data(simu.data)
#' n <- dim(simu.data)[1]
#' p <- dim(simu.data)[2]-2
#'
#' y <- simu.data[,1]
#' tr <- simu.data[,2]
#' x <- simu.data[,2+1:p]
#' x <- scale(x)
#'
#' # include only 10 covariates
#' x2 <- x[,1:10]
#'
#' ate.cal <- ate.nreg(y, tr, x2, ploss="cal", yloss="gaus")
#' matrix(unlist(ate.cal$est), ncol=2, byrow=TRUE,
#' dimnames=list(c("one", "ipw", "or", "est", "var", "ze", 
#' "diff.est", "diff.var", "diff.ze"), c("untreated", "treated")))
#' 
#' @export

ate.nreg <- function(y, tr, x, ploss="cal", yloss="gaus", off=NULL) {

  n <- length(tr)
  d <- 1+tr

  ps.nreg <- vector("list", 2)
  or.nreg <- vector("list", 2)

  mfp <- matrix(NA, nrow=n, ncol=2)
  mfo <- matrix(NA, nrow=n, ncol=2)

  for (k in 2:1) {

    # propensity score
    if (k==2 | ploss=="cal") {
      ps.nreg[[k]] <- glm.nreg(y= d==k, x=x, iw=NULL, loss=ploss)

      mfp[,k] <- ps.nreg[[k]] $fit

    } else { # k==1 & ploss=="ml"

      mfp[,1] <- 1-mfp[,2]
    }

    # outcome regression
    if (ploss=="ml") {
      iw <- NULL
    } else {  # "cal"
      iw <- 1/mfp[d==k,k] -1
    }

    or.nreg[[k]] <- glm.nreg(y[d==k], x[d==k,], iw=iw, loss=yloss)

    if (yloss=="gaus") {
      mfo[,k] <- c( cbind(1,x)%*%or.nreg[[k]] $coef )
    } else { # binary y
      mfo[,k] <- expit( c( cbind(1,x)%*%or.nreg[[k]] $coef ) )
    }
  }  

  # AIPW estimation
  out.aipw <- ate.aipw(y, tr, mfp=mfp, mfo=mfo, off=off)

  list(ps=ps.nreg, mfp=mfp, or=or.nreg, mfo=mfo, est=out.aipw)
}


###################################################################################
### non-regularized M-estimation using trust 
###################################################################################

# calibrated estimation
loss.cal <- function(gam, tr, x, iw)
{
  #gam: argument
  #tr: non-missingness/treatment indicator
  #x: covariate matrix; n*k
  
  n <- dim(x)[1]
  k <- dim(x)[2]
  
  if (!any(is.na(gam))) {
    w <- as.vector(x%*%gam)
    ew <- exp(-w)
    val <- mean(ifelse(tr, ew, w) *iw)
    
    grad <- apply(ifelse(tr, -ew, 1) *iw *x,2,mean)
    
    hess <- t(x)%*%(ifelse(tr,ew,0) *iw *x) /n
    
  } else {
    val <- Inf
    grad <- rep(NA, k)
    hess <- matrix(NA, k, k)
  }
  list(value=val, gradient=grad, hessian=hess)
}

# balance estimation
loss.bal <- function(gam, tr, x, iw)
{
  #gam: argument
  #tr: non-missingness/treatment indicator
  #x: constraint matrix; n*k
  
  n <- dim(x)[1]
  k <- dim(x)[2]
  
  if (!any(is.na(gam))) {
    g <- as.vector(x%*%gam)
    eg <- exp(-g)
    eg2 <- exp(g)

    val <- mean(ifelse(tr, eg-g, eg2+g) *iw)
    
    grad <- apply(ifelse(tr, -eg-1, eg2+1) *iw *x,2,mean)
    
    hess <- t(x)%*%(ifelse(tr,eg,eg2) *iw *x) /n
    
  }else {
    val <- Inf
    grad <- rep(NA, k)
    hess <- matrix(NA, k, k)
  }
  list(value=val, gradient=grad, hessian=hess)
}

#' Non-regularied M-estimation for fitting generalized linear models
#'
#' This function implements non-regularizd M-estimation for fitting generalized linear models with continuous or binary responses, including maximum likelihood, 
#' calibrated estimation, and covariate-balancing estimation in the latter case of fitting propensity score models. 
#'
#' Least squares estimation is implemented by calling \code{lm} for continuous responses (\code{loss}="gaus"). For binary responses, 
#' maximum likelihood estimation (\code{loss}="ml") is implemented by calling \code{glm}. Calibrated estimation (\code{loss}="cal") is implemented by 
#' using a trust-region algorithm in the R package \pkg{trust} to minimize the calibration loss, i.e., (8) in Tan (2017). 
#' Covariate-balancing estimation (\code{loss}="bal") in Imai and Ratkovic (2014) is implemented by using \pkg{trust} to minimize (38) in Tan (2017).
#'
#' @param y An \eqn{n} x \eqn{1} response vector.
#' @param x An \eqn{n} x \eqn{p} matix of covariates, excluding a constant.
#' @param iw An \eqn{n} x \eqn{1} weight vector. 
#' @param loss A loss function used, which can be specified as "gaus" for continuous responses, or "ml", "cal", or "bal" for binary responses.
#' @param init A \eqn{(p+1)} x \eqn{1} vector of initial values (the intercept and coefficients).
#'
#' @return
#' \item{coef}{The \eqn{(p+1)} x \eqn{1} vector of estimated intercept and coefficients.}
#' \item{fit}{The \eqn{n} x \eqn{1} vector of fitted values.}
#' \item{conv}{Logical; 1 if loss="gaus" for continuous responses or convergence is obtained within 1000 iterations by \code{glm} with loss="ml" or \code{trust} with loss="cal" or "bal" for binary responses.}
#'
#' @references
#' Imai, K. and Ratkovic, M. (2014) Covariate balancing propensity score, \emph{Journal of the Royal Statistical Society}, Ser. B, 76, 243-263.
#'
#' Tan, Z. (2017) Regularized calibrated estimation of propensity scores with model misspecification and high-dimensional data, arXiv:1710.08074. 
#' 
#' @examples
#' data(simu.data)
#' n <- dim(simu.data)[1]
#' p <- dim(simu.data)[2]-2
#'
#' y <- simu.data[,1]
#' tr <- simu.data[,2]
#' x <- simu.data[,2+1:p]
#' x <- scale(x)
#'
#' # include only 10 covariates
#' x2 <- x[,1:10]
#'
#' ps.ml <- glm.nreg(y=tr, x=x2, loss="ml")
#' check.ml <- mn.ipw(x2, tr, ps.ml$fit)
#' check.ml
#' 
#' ps.cal <- glm.nreg(y=tr, x=x2, loss="cal")
#' check.cal <- mn.ipw(x2, tr, ps.cal$fit)
#' check.cal  # should be numerically 0
#'
#' ps.bal <- glm.nreg(y=tr, x=x2, loss="bal")
#' check.bal <- mn.ipw(x2, tr, ps.bal$fit)
#' check.bal
#'
#' @export

glm.nreg <- function(y, x, iw=NULL, loss="cal", init=NULL) {

  if (is.null(iw)) {
    iw <- rep(1, length(y))
  } else {
    iw <- iw /mean(iw)
  }

  if (is.null(init)) 
    init <- lm(y~x, weights=iw)$coef

  x <- cbind(1,x)
   
  if (loss=="gaus") {
    solv <- lm(y ~ x[,-1, drop=F], weights=iw)
    coef <- solv$coef
    solv$conv <- 1

  } else if (loss=="ml") {
    solv <- glm(y ~ x[,-1, drop=F], family=binomial(link=logit), control=glm.control(maxit=1000), weights=iw)
    coef <- solv$coef

  } else if (loss=="cal") {
    solv <- trust(loss.cal, parinit=init, rinit=1, rmax=100, iterlim=1000, tr=y, x=x, iw=iw)
    coef <- solv$argument

  } else { #loss=="bal"
    solv <- trust(loss.bal, parinit=init, rinit=1, rmax=100, iterlim=1000, tr=y, x=x, iw=iw)
    coef <- solv$argument
  } 

  fit <- expit(as.vector(x%*%coef))
  
  list(coef=coef, fit=fit, conv=solv$conv)
}

###################################################################################
### regularized estimation with cross validation
###################################################################################

#' Regularied M-estimation for fitting generalized linear models based on cross validation
#'
#' This function implements regularized M-estimation for fitting generalized linear models with binary or contiunous responses 
#' based on cross validation.
#'
#' Cross validation is performed as described in Tan (2017, 2019). If not specified by users, the sequence of tuning parameters searched is defined as 
#' a geometric series of length \code{nrho}, starting from the value which yields a zero solution, and then decreasing by a factor \code{tune.fac} successively. 
#'
#' After cross validation, two tuning parameters are selected. The first and default choice is the value yielding the smallest average test loss.
#' The second choice is the largest value giving the average test loss within one standard error of the first choice (Hastie, Tibshirani, and Friedman 2016). 
#'
#' @param fold A fold number used for cross validation. 
#' @param nrho The number of tuning parameters searched in cross validation.
#' @param rho.seq A vector of tuning parameters searched in cross validation. If both \code{nrho} and \code{rho.seq} are specified, then \code{rho.seq} overrides \code{nrho}.
#' @param y An \eqn{n} x \eqn{1} response vector.
#' @param x An \eqn{n} x \eqn{p} matix of covariates, excluding a constant.
#' @param iw An \eqn{n} x \eqn{1} weight vector. 
#' @param loss A loss function, which can be specified as "guas" for continuous responses, or "ml" or "cal" for binary respones. 
#' @param n.iter The maximum number of iterations allowed as in \code{\link{glm.regu}}.
#' @param eps The tolerance used to declare convergence as in \code{\link{glm.regu}}. 
#' @param tune.fac The multiplier (factor) used to define \code{rho.seq} if only \code{nrho} is specified.
#' @param tune.cut Logical; if \code{TRUE}, all smaller tuning parameters are skipped once non-convergence is found with a tuning parameter. 
#' @param ann.init Logical; if \code{TRUE}, the estimates from the previous tuning parameter are used as the inital values when fitting with the current tuning parameter.
#' @param nz.lab A \eqn{p} x \eqn{1} logical vector (useful for simulations), indicating which covariates are included when calculating the number of nonzero coefficients. 
#' @param permut An \eqn{n} x \eqn{1} vector, giving a random permutation of the integers from 1 to \eqn{n}, which is used in cross validation. 
#'
#' @return
#' \item{permut}{An \eqn{n} x \eqn{1} vector, giving the random permutation used in cross validation.}
#' \item{rho}{The vector of tuning parameters, searched in cross validation.}
#' \item{non.conv}{A vector indicating the non-convergene status found or imputed if \code{tune.cut=TRUE}, for the tuning parmaters in cross validation.
#' For each tuning parameter, 0 indicates convergence, 1 non-convergence if exceeding \code{n.iter}, 2 non-convergence if exceeding \code{bt.lim}.}
#' \item{err.ave}{A vector giving the averages of the test losses in cross validation.}
#' \item{err.sd}{A vector giving the standard deviations of the test losses in cross validation.}
#' \item{sel.rho}{A vector of two selected tuning parameters by cross validation; see \strong{Details}.}
#' \item{sel.nz}{A vector of numbers of nonzero coefficients estimated for the selected tuning parameters.}
#' \item{sel.bet}{The \eqn{(p+1)} x \eqn{2} vector of estimated intercept and coefficients.}
#' \item{sel.fit}{The \eqn{n} x \eqn{2} vector of fitted values.}
#'
#' @references
#' Hastie, T., Tibshirani, R., and Friedman. J. (2016) \emph{The Elements of Statistical Learning} (second edition), Springer: New York.
#'
#' Tan, Z. (2017) Regularized calibrated estimation of propensity scores with model misspecification and high-dimensional data, arXiv:1710.08074. 
#'
#' Tan, Z. (2019) Model-assisted inference for treatment effects using regularized calibrated estimation with high-dimensional data, 
#' \emph{Annals of Statistics}, to appear (preprint arXiv:1801.09817).
#'
#' @examples 
#' \donttest{
#' data(simu.data)
#' n <- dim(simu.data)[1]
#' p <- dim(simu.data)[2]-2
#'
#' y <- simu.data[,1]
#' tr <- simu.data[,2]
#' x <- simu.data[,2+1:p]
#' x <- scale(x)
#'
#' ### Example 1: Regularized maximum likelihood estimation of propensity scores
#' ps.cv.rml <- glm.regu.cv(fold=5, nrho=1+10, y=tr, x=x, loss="ml")
#' ps.cv.rml$rho
#' ps.cv.rml$err.ave
#' ps.cv.rml$err.sd
#' ps.cv.rml$sel.rho
#' ps.cv.rml$sel.nz
#' 
#' fp.cv.rml <- ps.cv.rml $sel.fit[,1]
#' check.cv.rml <- mn.ipw(x, tr, fp.cv.rml)
#' check.cv.rml$est
#'
#' ### Example 2: Regularized calibrated estimation of propensity scores
#' ps.cv.rcal <- glm.regu.cv(fold=5, nrho=1+10, y=tr, x=x, loss="cal")
#' ps.cv.rcal$rho
#' ps.cv.rcal$err.ave
#' ps.cv.rcal$err.sd
#' ps.cv.rcal$sel.rho
#' ps.cv.rcal$sel.nz
#' 
#' fp.cv.rcal <- ps.cv.rcal $sel.fit[,1]
#' 
#' check.cv.rcal <- mn.ipw(x, tr, fp.cv.rcal)
#' check.cv.rcal$est
#' }
#' 
#' @export

glm.regu.cv <- function(fold, nrho=NULL, rho.seq=NULL, y, x, iw=NULL, loss="cal", 
                        n.iter=100, eps=1e-6, 
                        tune.fac=.5, tune.cut=TRUE, ann.init=TRUE,
                        nz.lab=NULL, permut=NULL) {

  n <- dim(x)[1]
  p <- dim(x)[2]

  cyc.seq <- NULL

  if (is.null(iw)) {
    iw <- rep(1, n)
  } else {
    iw <- iw /mean(iw)   ##<<-- normalized to have mean 1
  }

  # auto set rho.seq
  ybar <- mean(iw *y)

  if (loss=="gaus" || loss=="ml") {
    rho.max <- abs( apply(iw *(y-ybar) *x, 2, mean) )
  } else {
    rho.max <- abs( apply(iw *(y/ybar-1) *x, 2, mean) )
  }

  cyc.seq <- order(rho.max, decreasing=TRUE)

  if (is.null(rho.seq)) {
    if (is.null(nrho)) 
      stop("Error: nrho needs to be specified when rho.seq is not!")
 
    rho.max0 <- rho.max[cyc.seq[1]]
    rho.max0 <- my.signif(rho.max0, digits=3)

    rho.seq <- tune.fac^seq(nrho-1 ,0,-1) * rho.max0
    rho.seq <- signif(rho.seq, digits=3)
  }

  nrho <- length(rho.seq)  #rho.seq overiding nrho if both nrho and rho.seq are specified 

  #cat("rho.seq:", rho.seq, "\n")
  #cat("cyc.seq:", cyc.seq, "\n")

  #
  if (is.null(permut)) 
    permut <- sample(1:n)
  n.test <- n/fold

  err <- matrix(NA, nrow=nrho, ncol=fold)
  iter <- matrix(NA, nrow=nrho, ncol=fold)
  nz <- matrix(NA, nrow=nrho, ncol=fold)

  non.conv <- rep(0, nrho)

  #########################################
  for (j in 1:fold) {

    id.test <- (j-1)*n.test + 1:n.test
    test <- permut[id.test]
    train <- (1:n)[-test]

    wi.fits.i1 <- NULL

    for (i1 in nrho:1) {

      if (!non.conv[i1]) {

        rhos <- rep(rho.seq[i1], p)
 
        if (is.null(wi.fits.i1))  {
          out.regu <- glm.regu(y, x, iw, loss, init=NULL, rhos, test, offs=NULL, 
                               n.iter=n.iter, eps=eps, nz.lab=nz.lab)
        } else if (!ann.init) {
          out.regu <- glm.regu(y, x, iw, loss, init=NULL, rhos, test, offs=NULL,
                               wi.fits.i1$id, wi.fits.i1$Wmat, wi.fits.i1$Rmat, wi.fits.i1$zzs, wi.fits.i1$xxs,
                               n.iter=n.iter, eps=eps, nz.lab=nz.lab)
        } else {
          out.regu <- glm.regu(y, x, iw, loss, init=c(wi.fits.i1$inter, wi.fits.i1$bet), rhos, test, offs=NULL,
                               wi.fits.i1$id, wi.fits.i1$Wmat, wi.fits.i1$Rmat, wi.fits.i1$zzs, wi.fits.i1$xxs,
                               n.iter=n.iter, eps=eps, nz.lab=nz.lab)
        }

        # detect non-convergence 1=n.iter, 2=bt.lim
        if (out.regu$conv <=0) {
          if (tune.cut==1) { 
            non.conv[1:i1] <- 1 -out.regu$conv
          } else {
            non.conv[i1] <- 1 -out.regu$conv
          }

        } else {
          # use previous results
          wi.fits.i1 <- out.regu

          err[i1, j] <- out.regu$obj.test
          iter[i1, j] <- out.regu$iter
          nz[i1, j] <- out.regu$nz
        }
      } # if (!non.conv...)
    }  # for (i1 ...)

  }  # for (j ...)
  #########################################

  err.ave <- apply(err, 1, mean)
  err.sd <- apply(err, 1, sd) /sqrt(fold)

  iter.ave <- apply(iter, 1, mean)
  iter.max <- apply(iter, 1, max)

  nz.ave <- apply(nz, 1, mean)
  nz.max <- apply(nz, 1, max)

  # select the smallest error
  sel <- which.min(err.ave)
  sel <- c(sel, max(which(err.ave<= err.ave[sel]+err.sd[sel])))
  n.sel <- 2

  #
  iter.all <- rep(NA, nrho)
  nz.all <- rep(NA, nrho)

  sel.bet.all <- matrix(NA,1+p,n.sel)
  #sel.bet.nz <- vector("list", n.sel)
  sel.fit.all <- matrix(NA,n,n.sel)

  wi.fits.i1 <- NULL

  #########################################
  for (i1 in nrho:1) {

    sel.lab <- i1-sel==0
 
    if (any(sel.lab)) {

      rhos <- rep(rho.seq[i1], p)
           
      if (is.null(wi.fits.i1))  {
        out.regu <- glm.regu(y, x, iw, loss, init=NULL, rhos, test=NULL, offs=NULL, 
                             n.iter=n.iter, eps=eps, nz.lab=nz.lab)
      } else if (!ann.init) {
        out.regu <- glm.regu(y, x, iw, loss, init=NULL, rhos, test=NULL, offs=NULL, 
                             wi.fits.i1$id, wi.fits.i1$Wmat, wi.fits.i1$Rmat, wi.fits.i1$zzs, wi.fits.i1$xxs,
                             n.iter=n.iter, eps=eps, nz.lab=nz.lab)
      } else {
        out.regu <- glm.regu(y, x, iw, loss, init=c(wi.fits.i1$inter,wi.fits.i1$bet), rhos, test=NULL, offs=NULL, 
                             wi.fits.i1$id, wi.fits.i1$Wmat, wi.fits.i1$Rmat, wi.fits.i1$zzs, wi.fits.i1$xxs,
                             n.iter=n.iter, eps=eps, nz.lab=nz.lab)
      }

      if (out.regu$conv <=0) {
        stop("Error: conv<=0 for selected rho with the entire data")

      } else {
        # use previous results
        wi.fits.i1 <- out.regu

        iter.all[i1] <- out.regu$iter
        nz.all[i1] <- out.regu$nz

        sel.bet.all[, sel.lab] <- c(out.regu$inter, out.regu$bet)
        #print(sel.bet.all[, sel.lab, drop=F])
        sel.fit.all[, sel.lab] <- out.regu$fit

        #loc <- c(0, (1:p)[abs(out.regu$bet)>0])
        #for (jj in (1:n.sel)[sel.lab])
        #  sel.bet.nz[[jj]] <- cbind(loc, c(out.regu$inter, out.regu$bet[loc]))
        #print(sel.bet.nz[sel.lab])
      }
    } # if (any(sel.lab))

  } # for (i1 ...)
  #########################################

  list(permut=permut, rho=rho.seq, non.conv=non.conv,
       err.ave=err.ave, err.sd=err.sd, 

       sel.rho=rho.seq[sel],
       sel.nz=nz.all[sel], 
       sel.bet=sel.bet.all,
       sel.fit=sel.fit.all)
}

###################################################################################
### regularized estimation along a regularization path
###################################################################################

#' Regularied M-estimation for fitting generalized linear models along a regularization path
#'
#' This function implements regularized M-estimation for fitting generalized linear models with binary or contiunous responses 
#' along a regularization path.
#'
#' If not specified by users, the sequence of tuning parameters (i.e., the regularization path) is defined as a geometric series 
#' of length \code{nrho}, starting from the value which yields a zero solution, and then decreasing by a factor \code{tune.fac} successively. 
#'
#' @param nrho The number of tuning parameters in a regularization path.
#' @param rho.seq A vector of tuning parameters in a regularization path. If both \code{nrho} and \code{rho.seq} are specified, then \code{rho.seq} overrides \code{nrho}.
#' @param y An \eqn{n} x \eqn{1} response vector.
#' @param x An \eqn{n} x \eqn{p} matix of covariates, excluding a constant.
#' @param loss A loss function, which can be specified as "guas" for continuous responses, or "ml" or "cal" for binary respones. 
#' @param iw An \eqn{n} x \eqn{1} weight vector. 
#' @param n.iter The maximum number of iterations allowed as in \code{\link{glm.regu}}.
#' @param eps The tolerance used to declare convergence as in \code{\link{glm.regu}}. 
#' @param tune.fac The multiplier (factor) used to define rho.seq if only \code{nrho} is specified.
#' @param tune.cut Logical; if \code{TRUE}, all smaller tuning parameters are skipped once non-convergence is found with a tuning parameter. 
#' @param ann.init Logical; if \code{TRUE}, the estimates from the previous tuning parameter are used as the inital value when fitting with the current tuning parameter.
#' @param nz.lab A \eqn{p} x \eqn{1} logical vector (useful for simulations), indicating which covariates are included when calculating the number of nonzero coefficients. 
#'
#' @return
#' \item{rho}{The vector of tuning parameters included in the regularization path.}
#' \item{non.conv}{A vector indicating the non-convergene status found or imputed if \code{tune.cut=TRUE}, along the regularization path. 
#' For each tuning parameter, 0 indicates convergence, 1 non-convergence if exceeding \code{n.iter}, 2 non-convergence if exceeding \code{bt.lim}.}
#' \item{nz.all}{A vector giving the numbers of nonzero coefficients estimated, along the regularization path.}
#' \item{bet.all}{A matrix giving estimated intercept and coefficients, column by column, along the regularization path.}
#' \item{fit.all}{A matrix giving fitted values, column by column, along the regularization path.}
#'
#' @references
#' Tan, Z. (2017) Regularized calibrated estimation of propensity scores with model misspecification and high-dimensional data, arXiv:1710.08074. 
#'
#' Tan, Z. (2019) Model-assisted inference for treatment effects using regularized calibrated estimation with high-dimensional data, 
#' \emph{Annals of Statistics}, to appear (preprint arXiv:1801.09817).
#'
#' @examples 
#' \donttest{
#' data(simu.data)
#' n <- dim(simu.data)[1]
#' p <- dim(simu.data)[2]-2
#'
#' y <- simu.data[,1]
#' tr <- simu.data[,2]
#' x <- simu.data[,2+1:p]
#' x <- scale(x)
#'
#' ### Example 1: linear regression
#' out.rgaus.path <- glm.regu.path(rho.seq=c(.01, .02, .05, .1, .2, .5), y=y[tr==1], x=x[tr==1,], 
#'                                 loss="gaus")
#'
#' # the estimated intercept and coefficients; the first 10 are shown
#' out.rgaus.path$bet.all[1:10,]
#'
#' ### Example 2: logistic regression using likelihood loss
#' out.rml.path <- glm.regu.path(rho.seq=c(.002, .005, .01, .02, .05, .1), y=tr, x=x, loss="ml")
#' out.rml.path$bet.all[1:10,]
#'
#' ### Example 3: logistic regression using calibration loss
#' out.rcal.path <- glm.regu.path(rho.seq=c(.005, .01, .02, .05, .1, .2), y=tr, x=x, loss="cal")
#' out.rcal.path$bet.all[1:10,]
#' }
#'
#' @export

glm.regu.path <- function(nrho=NULL, rho.seq=NULL, y, x, iw=NULL, loss="cal",
                          n.iter=100, eps=1e-6, 
                          tune.fac=.5, tune.cut=TRUE, ann.init=TRUE,
                          nz.lab=NULL) {

  n <- dim(x)[1]
  p <- dim(x)[2]

  cyc.seq <- NULL

  if (is.null(iw)) {
    iw <- rep(1, n)
  } else {
    iw <- iw /mean(iw)   ##<<-- normalized to have mean 1
  }

  # auto set rho.seq
  ybar <- mean(iw *y)

  if (loss=="gaus" || loss=="ml") {
    rho.max <- abs( apply(iw *(y-ybar) *x, 2, mean) )
  } else {
    rho.max <- abs( apply(iw *(y/ybar-1) *x, 2, mean) )
  }

  cyc.seq <- order(rho.max, decreasing=TRUE)

  if (is.null(rho.seq)) {
    if (is.null(nrho)) 
      stop("Error: nrho needs to be specified when rho.seq is not!")
 
    rho.max0 <- rho.max[cyc.seq[1]]
    rho.max0 <- my.signif(rho.max0, digits=3)

    rho.seq <- tune.fac^seq(nrho-1, 0,-1) * rho.max0
    rho.seq <- signif(rho.seq, digits=3)
  }

  nrho <- length(rho.seq)  #rho.seq overiding nrho if both nrho and rho.seq are specified 

  #cat("rho.seq:", rho.seq, "\n")
  #cat("cyc.seq:", cyc.seq, "\n")

  #
  iter.all <- rep(NA, nrho)
  nz.all <- rep(NA, nrho)

  bet.all <- matrix(NA,1+p,nrho)
  #bet.nz <- vector("list", nrho)
  fit.all <- matrix(NA,n,nrho)

  non.conv <- rep(0, nrho)

  wi.fits.i1 <- NULL

  #########################################
  for (i1 in nrho:1) {

    if (!non.conv[i1]) {

      rhos <- rep(rho.seq[i1], p)
         
      if (is.null(wi.fits.i1))  {
        out.regu <- glm.regu(y, x, iw, loss, init=NULL, rhos, test=NULL, offs=NULL, 
                             n.iter=n.iter, eps=eps, nz.lab=nz.lab)
      } else if (!ann.init) {
        out.regu <- glm.regu(y, x, iw, loss, init=NULL, rhos, test=NULL, offs=NULL, 
                             wi.fits.i1$id, wi.fits.i1$Wmat, wi.fits.i1$Rmat, wi.fits.i1$zzs, wi.fits.i1$xxs,
                             n.iter=n.iter, eps=eps, nz.lab=nz.lab)
      } else {
        out.regu <- glm.regu(y, x, iw, loss, init=c(wi.fits.i1$inter,wi.fits.i1$bet), rhos, test=NULL, offs=NULL, 
                             wi.fits.i1$id, wi.fits.i1$Wmat, wi.fits.i1$Rmat, wi.fits.i1$zzs, wi.fits.i1$xxs,
                             n.iter=n.iter, eps=eps, nz.lab=nz.lab)
      }

      # detect non-convergence 1=n.iter, 2=bt.lim
      if (out.regu$conv <=0) {
        if (tune.cut==1) { 
          non.conv[1:i1] <- 1 -out.regu$conv
        } else {
          non.conv[i1] <- 1 -out.regu$conv
        }

      } else {
        # use previous results
        wi.fits.i1 <- out.regu

        iter.all[i1] <- out.regu$iter
        nz.all[i1] <- out.regu$nz

        bet.all[, i1] <- c(out.regu$inter, out.regu$bet)
        #print(bet.all[,i1, drop=F])
        fit.all[, i1] <- out.regu$fit

        #loc <- c(0, (1:p)[abs(out.regu$bet)>0])
        #bet.nz[[i1]] <- cbind(loc, c(out.regu$inter, out.regu$bet[loc]))
        #print(bet.nz[[i1]])
      }
    }

  } # for (i1 ...)
  #########################################

  list(rho=rho.seq, non.conv=non.conv,
       nz.all=nz.all, 
       bet.all=bet.all,
       fit.all=fit.all)
}

###################################################################################
### active set descent algorithm for linear and logistic (ml and cal) Lasso
###################################################################################

#' Regularied M-estimation for fitting generalized linear models with a fixed tuning parameter
#'
#' This function implements regularized M-estimation for fitting generalized linear models with continuous or binary responses 
#' for a fixed choice of tuning parameters.
#'
#' For continuous responses, this function uses an active-set descent algorithm (Osborne et al. 2000; Yang and Tan 2018) to solve the least-squares Lasso problem. 
#' For binary responses, regularized calibrated estimation is implemented using the Fisher scoring descent algorithm in Tan (2017), whereas
#' regularized maximum likelihood estimation is implemented in a similar manner based on quadratic approximation as in the R package \pkg{glmnet}.
#'
#' @param y An \eqn{n} x \eqn{1} response vector.
#' @param x An \eqn{n} x \eqn{p} matix of covariates, excluding a constant.
#' @param iw An \eqn{n} x \eqn{1} weight vector. 
#' @param loss A loss function, which can be specified as "guas" for continuous responses, or "ml" or "cal" for binary respones. 
#' @param init A \eqn{(p+1)} x \eqn{1} vector of initial values (the intercept and coefficients).
#' @param rhos A \eqn{p} x \eqn{1} vector of Lasso tuning parameters, usually a constant vector, associated with the \eqn{p} coefficients. 
#' @param test A vector giving the indices of observations between 1 and \eqn{n} which are included in the test set. 
#' @param offs An \eqn{n} x \eqn{1} vector of offset values, similarly as in \code{glm}.
#' @param id An argument which can be used to speed up computation.
#' @param Wmat An argument which can be used to speed up computation.
#' @param Rmat An argument which can be used to speed up computation.
#' @param zzs An argument which can be used to speed up computation.
#' @param xxs An argument which can be used to speed up computation.
#' @param n.iter The maximum number of iterations allowed. An iteration is defined by computing an quadratic approximation and solving a least-squares Lasso problem. 
#' @param eps The tolerance at which the difference in the objective (loss plus penalty) values is considered close enough to 0 to declare convergence.
#' @param bt.lim The maximum number of backtracking steps allowed.
#' @param nz.lab A \eqn{p} x \eqn{1} logical vector (useful for simulations), indicating which covariates are included when calculating the number of nonzero coefficients. 
#' If \code{nz.lab=NULL}, then \code{nz.lab} is reset to a vector of 0s. 
#' @param pos A value which can be used to facilitate recording the numbers of nonzero coefficients with or without the restriction by \code{nz.lab}. 
#' If \code{nz.lab=NULL}, then \code{pos} is reset to 1.
#'
#' @return
#' \item{iter}{The number of iterations performed up to \code{n.iter}.}
#' \item{conv}{1 if convergence is obtained, 0 if exceeding the maximum number of iterations, or -1 if exceeding maximum number of backtracking steps.}
#' \item{nz}{A value defined as (nz0 * \code{pos} + nz1) to record the numbers of nonzero coefficients without or with the restriction 
#' (denoted as nz0 and nz1) by \code{nz.lab}.
#' If \code{nz.lab=NULL}, then nz1 is 0, \code{pos} is 1, and hence \code{nz} is nz0.}
#' \item{inter}{The estimated intercept.}
#' \item{bet}{The \eqn{p} x \eqn{1} vector of estimated coefficients, excluding the intercept.}
#' \item{fit}{The vector of fitted values in the training set.}
#' \item{eta}{The vector of linear predictors in the training set.}
#' \item{tau}{The \eqn{p} x \eqn{1} vector of generalized signs, which should be -1 or 1 for a negative or positive estimate and between -1 and 1 for a zero estimate.}
#' \item{obj.train}{The average loss in the training set.}
#' \item{pen}{The Lasso penalty of the estimates.}
#' \item{obj}{The average loss plus the Lasso penalty.}
#' \item{fit.test}{The vector of fitted values in the test set.}
#' \item{eta.test}{The vector of linear predictors in the test set.}
#' \item{obj.test}{The average loss in the test set.}
#' \item{id}{This can be re-used to speed up computation.}
#' \item{Wmat}{This can be re-used to speed up computation.}
#' \item{Rmat}{This can be re-used to speed up computation.}
#' \item{zzs}{This can be re-used to speed up computation.}
#' \item{xxs}{This can be re-used to speed up computation.}
#'
#' @references
#' Osborne, M., Presnell, B., and Turlach, B. (2000) A new approach to variable selection in least squares problems, \emph{IMA Journal of Numerical Analysis}, 20, 389-404.
#'
#' Tan, Z. (2017) Regularized calibrated estimation of propensity scores with model misspecification and high-dimensional data, arXiv:1710.08074. 
#'
#' Yang, T. and Tan, Z. (2018) Backfitting algorithms for total-variation and empirical-norm penalized additive modeling with high-dimensional data, \emph{Stat}, 7, e198.
#'
#' Tibshirani, R. (1996) Regression shrinkage and selection via the Lasso, \emph{Journal of the Royal Statistical Society}, Ser. B, 58, 267-288.
#'
#' @examples 
#' data(simu.data)
#' n <- dim(simu.data)[1]
#' p <- dim(simu.data)[2]-2
#'
#' y <- simu.data[,1]
#' tr <- simu.data[,2]
#' x <- simu.data[,2+1:p]
#' x <- scale(x)
#'
#' ### Example 1: linear regression
#' # rhos should be a vector of length p, even though a constant vector
#' out.rgaus <- glm.regu(y[tr==1], x[tr==1,], rhos=rep(.05,p), loss="gaus")
#' 
#' # the intercept
#' out.rgaus$inter
#'
#' # the estimated coefficients and generalized signs; the first 10 are shown
#' cbind(out.rgaus$bet, out.rgaus$tau)[1:10,]
#'
#' # the number of nonzero coefficients 
#' out.rgaus$nz
#'
#' ### Example 2: logistic regression using likelihood loss
#' out.rml <- glm.regu(tr, x, rhos=rep(.01,p), loss="ml")
#' out.rml$inter
#' cbind(out.rml$bet, out.rml$tau)[1:10,]
#' out.rml$nz
#'
#' ### Example 3: logistic regression using calibration loss
#' out.rcal <- glm.regu(tr, x, rhos=rep(.05,p), loss="cal")
#' out.rcal$inter
#' cbind(out.rcal$bet, out.rcal$tau)[1:10,]
#' out.rcal$nz
#'
#' @export

glm.regu <- function(y, x, iw=NULL, loss="cal", init=NULL, rhos, test=NULL, offs=NULL, 
                     id=NULL, Wmat=NULL, Rmat=NULL, zzs=NULL, xxs=NULL,
                     n.iter=100, eps=1e-6, bt.lim=3, 
                     nz.lab=NULL, pos=10000) {

  n <- dim(x)[1]
  p <- dim(x)[2]

  if (is.null(iw)) {
    iw <- rep(1, n)
  } else {
    iw <- iw /mean(iw)   ##<<-- normalized to have mean 1
  }

  if (is.null(offs)) 
    offs <- rep(0,n) 

  if (is.null(nz.lab)) {
    nz.lab <- rep(0, p)
    pos <- 1
  }

  #
  if(!is.null(test)) {
    y.test <- y[test] 
    x.test <- x[test,] 
    offs.test <- offs[test]
    iw.test <- iw[test]         ##<<-- not normalized 

    y <- y[-test]
    x <- x[-test,] 
    offs <- offs[-test]
    iw <- iw[-test] /mean(iw[-test])  ##<<-- normalized

    n <- dim(x)[1]
  }

  obj.test <- NULL
  fit.test <- NULL
  eta.test <- NULL

  #
  xbar <- apply(iw *x,2,mean)
  x.til <- x - rep(xbar, each=n)

  if (!is.null(test))
    x.til.test <- x.test - rep(xbar, each=length(test))

  #
 #########################################
 if (loss=="gaus") {

  if (is.null(init)) {
    bets <- rep(0,p)  
  } else {
    bets <- init[-1]
  }

  y.til <- y -offs

  ybar <- mean(iw *y.til)
  y.til <- y.til -ybar

  if (is.null(id)) {
     out.asd <- asd(y.til, x.til, iw, bets, rhos=rhos,
                    nz.lab=nz.lab)
  } else {
     out.asd <- asd(y.til, x.til, iw, bets, rhos=rhos,
                    id, Wmat, Rmat, zzs, yxs=NULL, xxs,
                    nz.lab=nz.lab)
  } 

  nz <- out.asd$nz
  nz2 <- out.asd$nz2

  bets <- out.asd$bet

  inter0 <- sum(xbar*bets)
  inter <- ybar -inter0        # intercept before centering x

  #
  eta.til <- c(x.til%*%bets)
  eta <- offs +ybar +eta.til
  fit <- eta

  obj.train <- mean(iw *(y-fit)^2)/2
  pen <- sum(rhos*abs(bets)) 
  obj <- obj.train +pen

  if (!is.null(test)) {
    eta.test <- offs.test +ybar +c(x.til.test%*%bets)
    fit.test <- eta.test
    obj.test <- mean(iw.test *(y.test-fit.test)^2)/2
  }

  #
  taus <- apply(iw *(y-fit)*x.til, 2, mean)

  list(iter=1, conv=1, nz=nz *pos +nz2,
       inter=inter, bet=bets, fit=fit, eta=eta, tau=ifelse(rhos>0, taus/rhos, taus), obj.train=obj.train, pen=pen, obj=obj,
       fit.test=fit.test, eta.test=eta.test, obj.test=obj.test, 
       id=out.asd$id, Wmat=out.asd$Wmat, Rmat=out.asd$Rmat, zzs=out.asd$zzs, xxs=out.asd$xxs)

 #########################################
 } else if (loss=="ml") {

  ybar <- mean(iw *y)

  if (is.null(init)) {
    bets <- rep(0,p)
    inter <- log(ybar/(1-ybar))
  } else {
    inter <- init[1]
    bets <- init[-1]
  }

  bt <- 0

  w0 <- 1/4

  #
  inter0 <- sum(xbar*bets)
  inter <- inter +inter0  # this inter gives intercept after centering x; compare below

  eta.til <- c(x.til%*%bets)
  eta <- offs +inter +eta.til
  fit <- expit(eta)

  obj.train <- mean(iw *(-log(fit)+eta*(1-y)))
  pen <- sum(rhos*abs(bets)) 
  obj <- obj.train +pen

  #
  conv <- 0
  i <- 2   # not 1

  while (i<=n.iter && conv <= 0 && bt <= bt.lim) { # note conv <= 0

    #cat("####### i=", i, "\n")

    bets.old <- bets
    eta.old <- eta
    fit.old <- fit

    if (conv >= 0)
      wt <- rep(w0, n)

    y.til <- eta +(y-fit)/wt -offs   # -offs needed
  
    ywbar <- mean(iw *y.til *wt) /mean(iw *wt)
    y.til <- y.til -ywbar  # not ybar

    #
    if (i==2) {
      if (is.null(id)) {
        out.asd <- asd(y.til, x.til, iw, bets, rhos=rhos /wt[1], 
                       nz.lab=nz.lab)
      } else {
        out.asd <- asd(y.til, x.til, iw, bets, rhos=rhos /wt[1], 
                       id, Wmat, Rmat, zzs, yxs=NULL, xxs,
                       nz.lab=nz.lab)
      } 
    } else {
      out.asd <- asd(y.til, x.til, iw, bets, rhos=rhos /wt[1], 
                     id=out.asd$id, out.asd$Wmat, out.asd$Rmat, out.asd$zzs, yxs=NULL, out.asd$xxs,  # no need to multiply wt
                     nz.lab=nz.lab)
    }

    #if (i %% 10 ==2) {
    #  cat("i=", i, "\n")
    #  print(out.asd$nz)
    #  print(out.asd$bet)
    #}

    nz <- out.asd$nz
    nz2 <- out.asd$nz2

    bets <- out.asd$bet

    inter0 <- sum(xbar*bets)
    inter <- ywbar -inter0   # ywbar, not ybar; this inter gives intercept before centering x; then ywbar = inter + inter0; compare above

    #
    eta.til <- c(x.til%*%bets)
    eta <- offs +ywbar +eta.til   # ywbar, not inter, not ybar
    fit <- expit(eta)

    obj.train <- mean(iw *(-log(fit)+eta*(1-y)))
    pen <- sum(rhos*abs(bets))
    obj <- c(obj, obj.train +pen)

    #
    taus <- apply(iw *(y-fit)*x.til, 2, mean)

    if (obj[i-1]-obj[i] < -1e-6) {
      #stop("Error: non-decreasing objective!")
      conv <- -1
      wt <- wt *2
      bt <- bt+1

      if (i+1 <=n.iter && bt <= bt.lim) {
        obj <- c(obj[-i], obj[i-1])

        bets <- bets.old
        eta <- eta.old
        fit <- fit.old
      }
    } else {
      conv <- 0

      if  (obj[i-1]-obj[i] < eps)
        conv <- 1
    }

    i <- i+1
  } # end: while()

  if (!is.null(test)) {
    eta.test <- offs.test +ywbar +c(x.til.test%*%bets)
    fit.test <- expit(eta.test)
    obj.test <- mean(iw.test *(-log(fit.test)+eta.test*(1-y.test)))
  }

  list(iter=i-1, conv=conv, nz=nz *pos +nz2,
       inter=inter, bet=bets, fit=fit, eta=eta, tau=ifelse(rhos>0, taus/rhos, taus), obj.train=obj.train, pen=pen, obj=obj[i-1],
       fit.test=fit.test, eta.test=eta.test, obj.test=obj.test,
       id=out.asd$id, Wmat=out.asd$Wmat, Rmat=out.asd$Rmat, zzs=out.asd$zzs, xxs=out.asd$xxs)

 #########################################
 } else {  #loss=="cal"

  ybar <- mean(iw *y)

  if (is.null(init)) {
    bets <- rep(0,p)
    inter <- log(ybar/(1-ybar))
  } else {
    inter <- init[1]
    bets <- init[-1]
  }

  bt <- 0
 
  w0 <- 1

  #
  inter0 <- sum(xbar*bets)
  inter <- inter +inter0  # this inter gives intercept after centering x

  eta.til <- c(x.til%*%bets)
  eta <- offs +inter +eta.til
  fit <- expit(eta)

  obj.train <- mean(iw *ifelse(y, exp(-eta), eta))
  pen <- sum(rhos*abs(bets))
  obj <- obj.train +pen

  #
  conv <- 0
  i <- 2   # not 1

  while (i<=n.iter && conv <= 0 && bt <= bt.lim) { # note conv <= 0

    #cat("####### i=", i, "\n")

    bets.old <- bets
    eta.old <- eta
    fit.old <- fit

    if (conv >= 0)
      wt <- rep(w0, n)

    y.til <- eta +ifelse(y, 1/fit-1, -1)/wt -offs 

    ywbar <- mean(iw *y.til *wt) /mean(iw *wt)
    y.til <- y.til -ywbar  # not ybar

    #
    if (i==2) {
      if (is.null(id)) { 
        out.asd <- asd(y.til, x.til, iw, bets, rhos=rhos /wt[1],
                       nz.lab=nz.lab)
      } else {
        out.asd <- asd(y.til, x.til, iw, bets, rhos=rhos /wt[1], 
                       id, Wmat, Rmat, zzs, yxs=NULL, xxs,
                       nz.lab=nz.lab)
       } 
    } else {
      out.asd <- asd(y.til, x.til, iw, bets, rhos=rhos /wt[1],
                     id=out.asd$id, out.asd$Wmat, out.asd$Rmat, out.asd$zzs, yxs=NULL, out.asd$xxs,  # no need to multiply wt
                     nz.lab=nz.lab)
    }
   
    #if (i %% 10 ==2) {
    #  cat("i=", i, "\n")
    #  print(out.asd$nz)
    #  print(out.asd$bet)
    #}

    nz <- out.asd$nz
    nz2 <- out.asd$nz2

    bets <- out.asd$bet

    inter0 <- sum(xbar*bets)
    inter <- ywbar -inter0   # ywbar, not ybar; this inter gives intercept before centering x; then ywbar = inter + inter0; compare above

    #
    eta.til <- c(x.til%*%bets)
    eta <- offs +ywbar +eta.til   # ywbar, not inter, not ybar
    fit <- expit(eta)

    obj.train <- mean(iw *ifelse(y, exp(-eta), eta))
    pen <- sum(rhos*abs(bets))
    obj <- c(obj, obj.train +pen)

    #
    eta.til <- c(x.til%*%bets)
    eta.til.norm <- sqrt(mean(iw *eta.til^2))

    taus <- apply(iw *(y/fit-1)*x.til, 2, mean)

    if (obj[i-1]-obj[i] < -1e-6) {
      #stop("Error: non-decreasing objective!")
      conv <- -1
      wt <- wt *2
      bt <- bt+1

      if (i+1 <=n.iter && bt <= bt.lim) {
        obj <- c(obj[-i], obj[i-1])

        bets <- bets.old
        eta <- eta.old
        fit <- fit.old
      }

    } else {
      conv <- 0

      if  (obj[i-1]-obj[i] < eps)
        conv <- 1
    }

    i <- i+1
  } # end: while()

  if (!is.null(test)) {
    eta.test <- offs.test +ywbar +c(x.til.test%*%bets)
    fit.test <- expit(eta.test)
    obj.test <- mean(iw.test *ifelse(y.test, exp(-eta.test), eta.test))
  }

  list(iter=i-1, conv=conv, nz=nz *pos +nz2,
       inter=inter, bet=bets, fit=fit, eta=eta, tau=ifelse(rhos>0, taus/rhos, taus), obj.train=obj.train, pen=pen, obj=obj[i-1], 
       fit.test=fit.test, eta.test=eta.test, obj.test=obj.test,
       id=out.asd$id, Wmat=out.asd$Wmat, Rmat=out.asd$Rmat, zzs=out.asd$zzs, xxs=out.asd$xxs)

 #########################################
 } # end: loss=="cal"
}


###################################################################################
### active set descent algorithm for linear Lasso
###################################################################################

shift.down <- function(k, st, nz, id, Wmat, Rmat, zzs, yxs, xxs) {

  id2 <- id
  Wmat2 <- Wmat
  Rmat2 <- Rmat
  zzs2 <- zzs

  id2[st:(nz-1)] <- id[1+st:(nz-1)]
  id2[nz] <- id[st]

  for (j in st:(nz-1))
    if (j>1)
      Wmat2[1:(j-1),j] <- Wmat[1:j,j+1][-st]

  if (st>1) 
    Wmat2[1:(st-1), nz] <- Wmat[1:(st-1), st]
  Wmat2[st:(nz-1), nz] <- Wmat[st, (st+1):nz]

  if (nz<k) {
    for (j in (nz+1):k) {
      Wmat2[st:(nz-1),j] <- Wmat[(st+1):nz,j]
      Wmat2[nz,j] <- Wmat[st,j]
    }
  }

  #
  for (j in st:k) {

    if (j==1) {

      zzs2[1] <- xxs[id2[1]]

    } else {

      xzs <-  c(t(Rmat2[1:(j-1),1:(j-1)])%*% Wmat2[1:(j-1),j])

      b <- - xzs/zzs2[1:(j-1)]
      Rmat2[1:(j-1),j] <- Rmat2[1:(j-1),1:(j-1)]%*%b

      zzs2[j] <- xxs[id2[j]]+sum(b *xzs)
    } 
  }

  list(id=id2, Wmat=Wmat2, Rmat=Rmat2, zzs=zzs2)
}

shift.up <- function(k, iz, en, id, Wmat, Rmat, zzs, yxs, xxs) {

  id2 <- id
  Wmat2 <- Wmat
  Rmat2 <- Rmat
  zzs2 <- zzs

  id2[1+iz:(en-1)] <- id[iz:(en-1)]
  id2[iz] <- id[en]

  for (j in iz:(en-1))
    if (j>1)
      Wmat2[1:j,j+1][-iz] <- Wmat[1:(j-1),j]

  if (iz>1) 
    Wmat2[1:(iz-1), iz] <- Wmat[1:(iz-1), en]
  Wmat2[iz, (iz+1):en] <- Wmat[iz:(en-1), en]

  if (en<k) {
    for (j in (en+1):k) {
      Wmat2[(iz+1):en,j] <- Wmat[iz:(en-1),j]
      Wmat2[iz,j] <- Wmat[en,j]
    }
  }

  #
  for (j in iz:k) {

    if (j==1) {

      zzs2[1] <- xxs[id2[1]]

    } else {

      xzs <-  c(t(Rmat2[1:(j-1),1:(j-1)])%*% Wmat2[1:(j-1),j])

      b <- - xzs/zzs2[1:(j-1)]
      Rmat2[1:(j-1),j] <- Rmat2[1:(j-1),1:(j-1)]%*%b

      zzs2[j] <- xxs[id2[j]]+sum(b *xzs)
    } 
  }

  list(id=id2, Wmat=Wmat2, Rmat=Rmat2, zzs=zzs2)
}

patch <- function(x, k, en, id, Wmat, Rmat, zzs, yxs, xxs, iw) {

  if (en >k+1) {
    id2 <- id  

    id[k+1] <- id2[en]
    id[(k+2):en] <- id2[(k+1):(en-1)]
  }  

  if (k==0) {
    Wmat <- cbind(1)
    Rmat <- cbind(1)

    zzs <- xxs[id[k+1]]

  } else {
    Wmat <- rbind(cbind(Wmat, apply(iw *x[,id[1:k], drop=F]*x[,id[k+1]], 2, mean)),
                  c(rep(0,k), 1))

    xzs <-  c(t(Rmat)%*% Wmat[1:k,k+1])

    b <- - xzs/zzs
    Rmat <- rbind(cbind(Rmat, Rmat%*%b),
                  c(rep(0, k), 1))

    zzs <- c(zzs, xxs[id[k+1]]+sum(b *xzs))
  }

  list(id=id, Wmat=Wmat, Rmat=Rmat, zzs=zzs)
}


# Active-set descent algorithm for solving the least-squares Lasso problem
#
# This function implements an active-set descenet algorithm for solving the least-squares Lasso problem with a fixed choice of tuning parameters.
#
# @param y An \eqn{n} x \eqn{1} response vector.
# @param x An \eqn{n} x \eqn{p} matix of covariates, \emph{with columns centered}.
# @param iw An \eqn{n} x \eqn{1} weight vector. 
# @param init An initial vector of coefficients, excluding the intercept.
# @param rhos A \eqn{p} x \eqn{1} vector of tuning parameters, usually a constant vector, associated with the \eqn{p} coefficients. 
# @param id An argument which can be used to speed up computation.
# @param Wmat An argument which can be used to speed up computation.
# @param Rmat An argument which can be used to speed up computation.
# @param zzs An argument which can be used to speed up computation.
# @param yxs An argument which can be used to speed up computation.
# @param xxs An argument which can be used to speed up computation.
# @param nz.lab A \eqn{p} x \eqn{1} logical vector (useful for simulations), indicating which covariates are included when calculating the number of nonzero coefficients. 
#
# @return
# \item{cnt}{The number of active-set descent steps performed.}
# \item{k}{The number of covariates for which the Cholesky decomposision of the Gram matrix is computed.}
# \item{bet}{The vector of estimated intercept and coefficients.}
# \item{fit}{The vector of fitted values.}
# \item{tau}{The \eqn{p} x \eqn{1} vector of generalized signs, which should be -1 or 1 for a negative or positive and between -1 and 1 for a zero estimate.}
# \item{obj}{The average loss plus the Lasso penalty.}
# \item{id}{This can be re-used to speed up computation.}
# \item{Wmat}{This can be re-used to speed up computation.}
# \item{Rmat}{This can be re-used to speed up computation.}
# \item{zzs}{This can be re-used to speed up computation.}
# \item{yxs}{This can be re-used to speed up computation.}
# \item{xxs}{This can be re-used to speed up computation.}
# \item{nz2}{A value recording the numbers of nonzero coefficients with or without the restriction by \code{nz.lab}.}
#
# @references
# Osborne, M., Presnell, B., and Turlach, B. (2000) A new approach to variable selection in least squares problems, \emph{IMA Journal of Numerical Analysis}, 20, 389-404.
#
# Yang, T. and Tan, Z. (2018) Backfitting algorithms for total-variation and empirical-norm penalized additive modeling with high-dimensional data, \emph{Stat}, 7, e198.
# 
# @export

asd <- function(y, x, iw=NULL, init=NULL, rhos, 
                id=NULL, Wmat=NULL, Rmat=NULL, zzs=NULL, yxs=NULL, xxs=NULL, 
                nz.lab=NULL) {

  n <- dim(x)[1]
  p <- dim(x)[2]

  if (is.null(iw)) {
    iw <- rep(1, n)
  } else {
    iw <- iw /mean(iw)    ##<<-- normalized to have mean 1
  }

  if (is.null(nz.lab))
    nz.lab <- rep(0, p)

  #
  y.norm <- sqrt(mean(iw *y^2))

  if (is.null(yxs))
    yxs <- apply(iw *y*x, 2, mean)
  if (is.null(xxs)) 
    xxs <- apply(iw *x^2, 2, mean)
 
  #
  if (is.null(init)) {
    bets <- rep(0,p)
  } else {
    bets <- init
  }

  zgns <- sign(bets)
  nz <- sum(abs(zgns))

  #
  if (is.null(id)) {
    if (nz ==0) {
      id <- 1:p
    } else {
      id[1:nz] <- (1:p)[zgns!=0]
      if (nz<p)
        id[(nz+1):p] <- (1:p)[zgns==0]
    }
  }

  if (is.null(Wmat)) {
    k <- 0
  } else {
    k <- dim(Wmat)[1]
  }

  while (k < nz) {
    out.patch <- patch(x, k, k+1, id, Wmat, Rmat, zzs, yxs, xxs, iw)

    id <- out.patch$id
    Wmat <- out.patch$Wmat
    Rmat <- out.patch$Rmat
    zzs <- out.patch$zzs

    k <- k+1
  }

  #
  if (nz==0) {
    fit <- rep(0,n)
    fit.norm <- 0

    taus <- yxs
    obj <- y.norm^2/2

  } else {
    fit <- c(x[, id[1:nz], drop=F]%*%bets[id[1:nz]])
    fit.norm <- sqrt(mean(iw *fit^2))

    taus <- apply(iw *(y-fit)*x, 2, mean)
    obj <- mean(iw *(y-fit)^2)/2 +sum(rhos*abs(bets))
  }

  #
  flag <- 0
  cnt <- 0
  cnt2 <- c(0,0)

 if (all(abs(yxs)<=rhos)) {
  bets <- rep(0,p)
  zgns <- rep(0,p)
  nz <- 0

  fit <- rep(0,n)
  taus <- yxs
  obj <- c(obj, y.norm^2/2)

  flag <- -1

 } else {

  #########################################
  while (flag==0) {
    rho <- rhos[id]

    #
    z.id <- NULL
    del <- NULL

    if (nz>0) {
      bet <- bets[id] [1:nz]
      zgn <- zgns[id] [1:nz]

      # depend on zgn only 
      bet.til <- t(Rmat[1:nz,1:nz])%*%(yxs[id][1:nz]-rho[1:nz]*zgn) /ifelse(zzs[1:nz]> 1e-10, zzs[1:nz], 1e-10)
      bet2 <- c(Rmat[1:nz,1:nz]%*%bet.til)

      del <- 1
      for (j in 1:nz) {
        if (j==1) temp <- bet2[1]
        else temp <- (1-del)*bet[j]+del*bet2[j]

        if (rho[j]>0 && sign(temp)*zgn[j] <0) {
          z.id <- j
          del <- -bet[j] /(bet2[j]-bet[j])
        }
      }
    }
 
    #
    if (is.null(z.id)) {
  
      if (nz>0) {
        cnt <- cnt+1
        cnt2[1] <- cnt2[1]+1

        bets[id][1:nz] <- bet2

        # 
        fit <- c(x[, id[1:nz], drop=F]%*%bet2)
        fit.norm <- sqrt(mean(iw *fit^2))

        taus <- apply(iw *(y-fit)*x, 2, mean)
        obj <- c(obj, mean(iw *(y-fit)^2)/2 +sum(rhos*abs(bets)))
      }

      if (nz <p) {
        tau <- taus[id] [(nz+1):p]

          j <- which.max( (abs(tau)-rho[(nz+1):p])/sqrt(xxs[id][(nz+1):p]) )
  
          if (abs(tau[j]) -rho[nz+j] > 1e-10) {
            nz.id <- j
  
            # patch up
            if (nz+nz.id >k) {

              out.patch <- patch(x, k, nz+nz.id, id, Wmat, Rmat, zzs, yxs, xxs, iw)

              id <- out.patch$id
              Wmat <- out.patch$Wmat
              Rmat <- out.patch$Rmat
              zzs <- out.patch$zzs
 
              nz.id <- k+1 -nz 
              k <- k+1
            }

            zgns[id][nz+nz.id] <- sign(tau[j])

            if (nz.id>1) {
              out.shift <- shift.up(k, nz+1, nz+nz.id, id, Wmat, Rmat, zzs, yxs, xxs)

              id <- out.shift$id
              Wmat <- out.shift$Wmat
              Rmat <- out.shift$Rmat
              zzs <- out.shift$zzs
            }

            nz <- nz+1 

          } else {   # abs(tau[j])>rho[nz+j]
            flag <- 2
          }

      } else {   # nz<p
        flag <- 1
      }
    } else {   # is.null(z.id)

      cnt <- cnt+1
      cnt2[2] <- cnt2[2]+1

      bets[id][1:nz] <- bet + del*(bet2-bet)
      bets[id][z.id] <- 0
      zgns[id][z.id] <- 0

      if (z.id<nz) {
        out.shift <- shift.down(k, z.id, nz, id, Wmat, Rmat, zzs, yxs, xxs)

        id <- out.shift$id
        Wmat <- out.shift$Wmat
        Rmat <- out.shift$Rmat
        zzs <- out.shift$zzs
      }

      nz <- nz-1

      # 
      if (nz==0) {
        fit <- rep(0,n)
        fit.norm <- 0

        taus <- yxs
        obj <- y.norm^2/2
      } else {
        fit <- c(x[, id[1:nz], drop=F]%*%bets[id[1:nz]])
        fit.norm <- sqrt(mean(iw *fit^2))

        taus <- apply((y-fit)*x, 2, mean)
        obj <- c(obj, mean(iw *(y-fit)^2)/2 +sum(rhos*abs(bets)))
      }
    }

  #########################################
  }  # end: flag==0
 }  # end: all(abs(yxs)<=rhos)

  list(cnt=cnt, k=k,
       bet=bets, nz=nz, fit=fit, tau=ifelse(rhos>0, taus/rhos, taus), obj=obj, 
       id=id, Wmat=Wmat, Rmat=Rmat, zzs=zzs, yxs=yxs, xxs=xxs,
       nz2=sum(nz.lab[zgns!=0]))
}

# a documentation for the package

#' RCAL: Regularized calibrated estimation
#'
#' Regularized calibrated estimation for causal inference and missing-data problems with high-dimensional data.
#'
#' The R package \code{RCAL} - version 1.0 can be used for two main tasks:
#' \itemize{
#'   \item to estimate the mean of an outcome in the presence of missing data,
#'   \item to estimate the average treatment effects in causal inference.
#' }
#' There are 3 high-level functions provided for the first task:
#' \itemize{
#'   \item \code{mn.nreg}: inference using non-regularized calibrated estimation,
#'   \item \code{mn.regu.cv}: inference using regularized calibrated estimation based on cross validation,
#'   \item \code{mn.regu.path}: inference using regularized calibrated estimation along a regularization path. 
#' }
#' The first function \code{mn.nreg} is appropriate only in relatively low-dimensional settings,
#' whereas the functions \code{mn.regu.cv} and \code{mn.regu.path} are designed to deal with high-dimensional data (namely, 
#' the number of covariates close to or greater than the sample size). 
#' In parallel, there are 3 functions for the second task, \code{ate.nreg}, \code{ate.regu.cv}, and \code{ate.regu.path}.
#' These functions can also be used to perform inference for the average treatment effects on the treated or on the untreated.
#' Currently, the treatment is assumed to be binary (i.e., untreated or treated). 
#' Extensions to multi-valued treatments will be incorporated in later versions. 
#'
#' The package also provides lower-level functions, including \code{glm.nreg} to implement non-regularized M-estimation and \code{glm.regu} to 
#' implement Lasso regularized M-estimation for fitting generalized linear models currently with continuous or binary outcomes. 
#' The latter function \code{glm.regu} uses an active-set descent algorithm, which enjoys a finite termination property 
#' for solving least-squares Lasso problems. 
#'
#' See the the vignette for more details.
#' @docType package
#' @name RCAL-package
#' @aliases RCAL
NULL
