#' @title Estimete the parameters with gmm after IV selecting
#' @description Hybrid gmm estimator after selecting IVs in the reduced form equation.
#' @param g A function of the form \eqn{g(\theta,x)} and which returns a \eqn{n \times q} matrix with typical element \eqn{g_i(\theta,x_t)} for \eqn{i=1,...q} and \eqn{t=1,...,n}. This matrix is then used to build the q sample moment conditions. It can also be a formula if the model is linear  (see details gmm).
#' @param x The design matrix, without an intercept.
#' @param z The instrument variables matrix.
#' @param max.degree The upper limit value of degree of B-splines when using BIC/AIC to choose the tuning parameters, default is BIC.
#' @param criterion The criterion by which to select the regularization parameter. One of "AIC", "BIC", "GCV", "AICc","EBIC", default is "BIC".
#' @param df.method How should effective model parameters be calculated? One of: "active", which counts the number of nonzero coefficients; or "default", which uses the calculated df returned by grpreg, default is "default".
#' @param penalty The penalty to be applied to the model. For group selection, one of grLasso, grMCP, or grSCAD. For bi-level selection, one of gel or cMCP, default is " grLasso".
#' @param endogenous.index Specify which variables in design matrix are endogenous variables, the variable corresponds to the value 1 is endogenous variables,  the  variable corresponds to the value 0 is exogenous variable, the default is all endogenous variables.
#' @param IV.intercept Intercept of instrument variables, default is “FALSE”.
#' @param family Either "gaussian" or "binomial", depending on the response.default is " gaussian ".
#' @param ... Arguments passed to gmm (such as type,kernel...,detail see gmm).
#' @details See naivereg and gmm.
#' @return An object of type \code{naive.gmm} which is a list with the following
#' components:
#' \item{degree}{Degree of B-splines.}
#' \item{criterion}{The criterion by which to select the regularization parameter. One of "AIC", "BIC", "GCV", "AICc","EBIC", default is "BIC".}
#' \item{ind}{The index of selected instrument variables.}
#' \item{ind.b}{The index of selected instrument variables after B-splines.}
#' \item{gmm}{Gmm object, detail see gmm.}
#' @author Qingliang Fan, KongYu He, Wei Zhong
#' @references Q. Fan and W. Zhong (2018), “Nonparametric Additive Instrumental Variable Estimator: A Group Shrinkage Estimation Perspective,” Journal of Business & Economic Statistics, doi: 10.1080/07350015.2016.1180991.
#' @references Caner, M. and Fan, Q. (2015), Hybrid GEL Estimators: Instrument Selection with Adaptive Lasso, Journal of Econometrics, Volume 187, 256–274.
#' @examples
#'# gmm estimation after IV selection
#'data("naivedata")
#'x=naivedata[,1]
#'y=naivedata[,2]
#'z=naivedata[,3:22]
#'naive.gmm(y~x+x^2,cbind(x,x^2),z)
#' @export
naive.gmm <- function(g,x,z,max.degree=10,criterion=c("BIC","AIC","GCV","AICc","EBIC"),
                           df.method=c("default","active"),penalty = c("grLasso", "grMCP", "grSCAD", "gel", "cMCP"),endogenous.index=c(),IV.intercept=FALSE,
                           family=c("gaussian", "binomial", "poisson"),...) {
  x=as.matrix(x)
  z=as.matrix(z)
  if(dim(as.matrix(z))[2]<dim(as.matrix(x))[2]){
    stop('The number of selected instruments is less than the number of endogenous variables, please find other instruments')
  }
  if(max.degree <= 3){
    stop("max.degree is at least 3")
  }
  if(missing(criterion)){
    criterion='BIC'
  }
  if(missing(df.method)){
    df.method='default'
  }
  if(missing(penalty)){
    penalty='grLasso'
  }
  if(missing(endogenous.index)){
    endogenous.index=rep(1,dim(x)[2])
  }
  if(missing(family)){
    family='gaussian'
  }

  z=cbind(z,x[,endogenous.index!=1])
  x=x[,endogenous.index]


  Zselect<-IVselect(z,x,criterion=criterion,df.method=df.method,penalty=penalty,IV.intercept=IV.intercept,family = family)
  IV <- Zselect$IVselect
  cnames <-paste("x",1:dim(IV)[2],sep="")
  colnames(IV)<-cnames

  naivelasso.gmm<-gmm(g = g, x = IV,...)

  list(degree=Zselect$degree,criterion=Zselect$criterion,ind=Zselect$ind,ind.b=Zselect$ind.b,gel=naivelasso.gmm)

}
