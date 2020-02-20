#' @title Logistic-regression Instrumental Variables Estimator
#' @description Binary endogenous variables are commonly encountered in program evaluations using observational data. This is a two-stage approach to estimate the dummy endogenous treatment effect using high-dimensional instrumental variables (IV). In the first stage, we use a penalized logistic reduced form model to accommodate both the binary nature of the endogenous treatment and the high-dimensionality of instrumental variables. In the second stage, we replace the original treatment variable by its estimated propensity score and run a least squares regression to obtain a penalized Logistic-regression Instrumental Variables Estimator (LIVE). If the structural equation model is also high-dimensional, one could use DS-LIVE in this package for selecting both the control variables and IVs.
#' @param y Response variable, an N x 1 vector.
#' @param x The design matrix, including endogenous variable, the value of endogenous variable is 0 or 1 (binary).
#' @param z The instrumental variables matrix.
#' @param penalty The penalty to be applied to the model. Either "SCAD" (the default), "MCP", or "lasso".
#' @param nfolds The response number of folds - default is 5. Although nfolds can be as large as the sample size (leave-one-out CV), it is not recommended for large datasets. Smallest value allowable is nfolds=3.
#' @param endogenous.index Specify which variables in design matrix are endogenous variables, the variable corresponds to the value 1 is endogenous variables,  the  variable corresponds to the value 0 is exogenous variable, the default is all endogenous variables.
#' @param gamma The tuning parameter of the MCP/SCAD penalty. Default is 3.7.
#' @param alpha Tuning parameter for the Mnet estimator which controls the relative contributions from the MCP/SCAD penalty and the ridge, or L2 penalty. alpha=1 is equivalent to MCP/SCAD penalty, while alpha=0 would be equivalent to ridge regression. However, alpha=0 is not supported; alpha may be arbitrarily small, but not exactly 0.
#' @param lambda.min The smallest value for lambda, as a fraction of lambda.max, default is 0.05.
#' @param nlambda The number of lambda values, default is 100.
#' @param ... other arguments.
#' @details This is a two stage estimation. In the first stage, a high-dimensional logistic reduced form model with penalty (such as SCAD, lasso, etc.) is used to approximate the optimal instrument. In the second stage, we replace the original treatment variable by its estimated propensity score and run a least squares regression to obtain the penalized Logistic-regression Instrumental Variables Estimator (LIVE). The large dimensional IV could be the original variables or the functional transformations such as series, B-spline functions, etc.
#' @return An object of type \code{LIVE} which is a list with the following
#' components:
#' \item{coefficients}{The coefficients of x.}
#' \item{lambda.min}{The value of lambda that gives minimum cvm.}
#' \item{ind}{The selected variables of z.}
#' \item{Xhat}{The xhat estimated by z.}
#' \item{IVnum}{The number of instrumented variables after filtering.}
#' \item{penalty}{Same as above.}
#' \item{alpha}{Same as above.}
#' @author Qingliang Fan, KongYu He, Wei Zhong
#' @references Wei Zhong, Wei Zhou, Qingliang Fan and Yang Gao (2020), “Dummy Endogenous Treatment Effect Estimation Using High-Dimensional Instrumental Variables”, working paper.
#' @examples
#'#Logistic-regression Instrumental Variables Estimator
#'data("LIVEdata")
#'y=LIVEdata[,1]
#'x=LIVEdata[,2]
#'z=LIVEdata[,3:52]
#'res = LIVE(y,x,z,penalty='SCAD',gamma = 3.7,alpha = 1,lambda.min = 0.05)
#' @export


LIVE <- function(y,x,z,penalty=c("SCAD","MCP","lasso"),nfolds=5,endogenous.index=c(),
                 gamma = 3.7,alpha = 1,lambda.min = 0.05,nlambda = 100,...) {
  if(dim(as.matrix(z))[2]<dim(as.matrix(x))[2]){
    stop('The number of selected instruments is less than the number of endogenous variables, please find other instruments')
  }
  if(missing(endogenous.index)){
    endogenous.index=rep(1,dim(as.matrix(x))[2])
  }
  if(missing(penalty)){
    penalty='SCAD'
  }

  x=as.matrix(x)
  z=as.matrix(z)
  y=as.matrix(y)

  z=cbind(z,x[,endogenous.index!=1])
  x=x[,endogenous.index]

  cv.scad<-cv.ncvreg(z,x,nflods=nfolds)
  bestlam.scad<-cv.scad$lambda.min

  fit.hd<-ncvreg(z,x,family = 'binomial',penalty = penalty,gamma = gamma,alpha = alpha,lambda.min = lambda.min,nlambda = nlambda,lambda=bestlam.scad,eps = 1e-4,max.iter = 50000)
  coefscad<-fit.hd$beta
  coef.scad<-coefscad[-1]
  #which(coef.scad!=0)

  ind = which(coef.scad!=0)
  if(all(coefscad==0)){
    beta.live_scad=1
  }else{
    se.va.scad<-z[,ind]
    num.s<-ncol(se.va.scad)
    glm.scad=glm(x~se.va.scad[,1:num.s],family = 'binomial',control=list(maxit=1000))
    xhat.scad=glm.scad$fitted
    beta.live_scad=lm(y~xhat.scad-1)$coef
  }

  list(coefficients=beta.live_scad,lambda.min=bestlam.scad,ind=ind,Xhat=xhat.scad,IVnum = num.s,penalty = penalty,alpha = alpha)
}
