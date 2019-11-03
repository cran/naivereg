#' @title Logistic-regression Instrumental Variables Estimator
#' @description A two-stage approach to estimate endogenous treatment effects using high-dimensional instrumental variables.
#' @param y Response variable, an Nx1 vector.
#' @param x The design matrix, including endogenous variable, the value of endogenous variable is 0 or 1 (binary).
#' @param z The instrumental variables matrix.
#' @param family Quantitative for family="gaussian", or family="poisson" (non-negative counts). For family="binomial" should be either a factor with two levels, or a two-column matrix of counts or proportions (the second column is treated as the target class; for a factor, the last level in alphabetical order is the target class). For family="multinomial", can be a nc>=2 level factor, or a matrix with nc columns of counts or proportions. For either "binomial" or "multinomial", if y is presented as a vector, it will be coerced into a factor. For family="cox", y should be a two-column matrix with columns named 'time' and 'status'. The latter is a binary variable, with '1' indicating death, and '0' indicating right censored. The function Surv() in package survival produces such a matrix. For family="mgaussian", y is a matrix of quantitative responses.
#' @param penalty The penalty to be applied to the model. Either "SCAD" (the default), "MCP", or "lasso".
#' @param nfolds The response number of folds - default is 5. Although nfolds can be as large as the sample size (leave-one-out CV), it is not recommended for large datasets. Smallest value allowable is nfolds=3.
#' @param endogenous.index Specify which variables in design matrix are endogenous variables, the variable corresponds to the value 1 is endogenous variables,  the  variable corresponds to the value 0 is exogenous variable, the default is all endogenous variables.
#' @param gamma The tuning parameter of the MCP/SCAD penalty. Default is 3.7.
#' @param alpha Tuning parameter for the Mnet estimator which controls the relative contributions from the MCP/SCAD penalty and the ridge, or L2 penalty. alpha=1 is equivalent to MCP/SCAD penalty, while alpha=0 would be equivalent to ridge regression. However, alpha=0 is not supported; alpha may be arbitrarily small, but not exactly 0.
#' @param lambda.min The smallest value for lambda, as a fraction of lambda.max. Default is 0.05
#' @param nlambda The number of lambda values. Default is 100.
#' @param ... other arguments.
#' @details This is a two stage estimation. In the first stage, a high-dimensional logistic reduced form model with penalty (such as SCAD, lasso, etc.) is used to approximate the optimal instrument. In the second stage, we replace the original treatment variable by its estimated propensity score and run a least squares regression to obtain the penalized Logistic-regression Instrumental Variables Estimator (LIVE). The large dimensional IV could be the original variables or the functional transformations such as series, B-spline functions, etc.
#' @return An object of type \code{LIVE} which is a list with the following
#' components:
#' \item{coefficients}{The coefficients of x.}
#' \item{lambda.min}{The value of lambda that gives minimum cvm.}
#' \item{ind}{The selected variables of z.}
#' \item{Xhat}{The xhat estimated by z.}
#' \item{IVnum}{The number of instrumented variables after filtering.}
#' \item{family}{Same as above.}
#' \item{penalty}{Same as above.}
#' \item{alpha}{Same as above.}
#' @author Qingliang Fan, KongYu He, Wei Zhong
#' @references Wei Zhong, Wei Zhou, Qingliang Fan and Yang Gao (2019), “Estimating Endogenous Treatment Effect Using High-Dimensional Instruments with an Application to the Olympic Effect”
#' @examples
#'#Logistic-regression Instrumental Variables Estimator
#'data("LIVEdata")
#'y=LIVEdata[,1]
#'x=LIVEdata[,2]
#'z=LIVEdata[,3:52]
#'res = LIVE(y,x,z,family='binomial',penalty='SCAD',gamma = 3.7,alpha = 1,lambda.min = 0.05)
#' @export


LIVE <- function(y,x,z,family=c("gaussian","poisson","binomial","multinomial","cox","mgaussian"),penalty=c("SCAD","MCP","lasso"),nfolds=5,endogenous.index=c(),
                 gamma = 3.7,alpha = 1,lambda.min = 0.05,nlambda = 100,...) {
  if(dim(as.matrix(z))[2]<dim(as.matrix(x))[2]){
    stop('The number of selected instruments is less than the number of endogenous variables, please find other instruments')
  }
  if(missing(endogenous.index)){
    endogenous.index=rep(1,dim(as.matrix(x))[2])
  }
  if(missing(family)){
    stop("family is missing")
  }
  if(missing(penalty)){
    criterion='SCAD'
  }
  if(missing(family)){
    criterion='binomial'
  }
  x=as.matrix(x)
  z=as.matrix(z)
  Y=as.matrix(y)

  Z=cbind(z,x[,endogenous.index!=1])
  X=x[,endogenous.index]

  cv.scad<-cv.ncvreg(z,x,nflods=nfolds)
  bestlam.scad<-cv.scad$lambda.min

  fit.hd<-ncvreg(z,x,family = family,penalty = penalty,gamma = gamma,alpha = alpha,lambda.min = lambda.min,nlambda = nlambda,lambda=bestlam.scad,eps = 1e-4,max.iter = 10000)
  coefscad<-fit.hd$beta
  coef.scad<-coefscad[-1]
  #which(coef.scad!=0)

  ind = which(coef.scad!=0)
  if(all(coefscad==0)){
    beta.live_scad=1
  }else{
    se.va.scad<-z[,ind]
    num.s<-ncol(se.va.scad)
    glm.scad=glm(x~se.va.scad[,1:num.s],family = family)
    xhat.scad=glm.scad$fitted
    beta.live_scad=lm(y~xhat.scad-1)$coef
  }

  list(coefficients=beta.live_scad,lambda.min=bestlam.scad,ind=ind,Xhat=xhat.scad,IVnum = num.s,family = family,penalty = penalty,alpha = alpha)
}
