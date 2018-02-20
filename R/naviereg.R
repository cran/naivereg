#' @title Nonparametric additive instrumental variable estimator
#' @description NAIVE is the nonparametric additive instrumental variable estimator with the
#' adaptive group Lasso. It is using group lasso and B-splines to obtain the valid instrument
#' variables where BIC or EBIC are applied to choose the tuning parameters. Then we get the
#' two-stage least squares (2SLS) estimator with selected IV.
#' @param y Response variable, a matrix Nx1
#' @param x The design matrix, without an intercept
#' @param z The instrument variables matrix
#' @param max.degree The upper limit value of degree of B-splines when using BIC/AIC to choose the tuning parameters, default is BIC.
#' @param intercept Estimate with intercept or not, default is "TRUE"
#' @param criterion The criterion by which to select the regularization parameter. One of "AIC", "BIC", "GCV", "AICc", or "EBIC"; default is "BIC".
#' @param df.method How should effective model parameters be calculated? One of: "active", which counts the number of nonzero coefficients; or "default", which uses the calculated df returned by grpreg. default is "default".
#' @param penalty The penalty to be applied to the model. For group selection, one of grLasso, grMCP, or grSCAD. For bi-level selection, one of gel or cMCP. Default is " grLasso".
#' @param endogenous.index Specify which variables in design matrix are endogenous variables, the  variable corresponds to the value 1 is endogenous variables,  the  variable corresponds to the value 0  is exogenous variable, the default is all endogenous variables
#' @param IV.intercept Intercept of instrument variables, default is “FALSE”
#' @param family Either "gaussian" or "binomial", depending on the response.default is " gaussian "
#' @details Consider the following structural equation with endogenous regressors\eqn{Y_{i}= x_{u}^{T}\beta+ \epsilon_{i}}
#'
#'To solve the endogeneity problem, instrumental variables are employed to obtain a consistent estimator of the population regression coefficient \eqn{\beta}. In practice, many potential instruments, including their series terms, may be recruited to approximate the optimal instrument and improve the precision of IV estimators. On the other hand, if many irrelevant instruments are contained in the reduced form equation, the approximation of the optimal instrument is generally unsatisfactory and the IV estimator is less efficient. In some cases where the dimensionality of \eqn{z_{i}}is even higher than the sample size, the linear IV method fails. To address these issues, the model sparsity is usually assumed and the penalized approaches can be applied to improve the efficiency of IV estimators. In this paper we propose the first-stage parsimonious predictive models and estimate optimal instruments in IV models with potentially more instruments than the sample size n.
#'
#'The performance of the linear IV estimator in the finite sample is largely dependent on the validity of linearity assumption. This phenomenon motivated us to consider a more general nonlinear reduced form equation to capture as much information of \eqn{x_{i}}as possible using instruments \eqn{z_{i}} under the high-dimensional model settings. This nonparametric idea for the reduced form model is consistent with Newey (1990). We consider the following nonparametric additive reduced form model with a large number of possible instruments.
#'
#'\eqn{x_{il} = \mu_l+\sum_{j=1}^p f_{ij}z_{ij}+\xi_{il}}
#'
#'To estimate the nonparametric components above, we use B-spline basis functions by following the idea of Huang, Horowitz, and Wei (2010). Let \eqn{S_{n}}be the space of polynomial splines of degrees L>1 and let \eqn{\phi_{k},k=1,2,…,m_{n}}be normalized B-spline basis functions for \eqn{S_{n}}, where  \eqn{m_{n}} is the sum of the polynomial degree L and the number of knots. Let be the \eqn{\psi_{k}(z_{ij})=\phi_{k}(z_{ij})-n^{-1}\sum_{i=1}^n \phi_{k}(z_{ij})}centered B-spline basis functions for the  th instrument. The model can then be rewritten using an approximate linear reduced form:
#'
#'\eqn{x_{il} = \mu_{l}+\sum_{j=1}^pf_{ij} \sum_{k=1}^{m_{n}} (\gamma_{ij})\psi(z_{ij})+\xi_{il}}
#'
#'To select the significant instruments and estimate the component functions simultaneously, we consider the following penalized objective function with an adaptive group Lasso penalty (Huang, Horowitz, and Wei 2010) for each \eqn{l} th endogenous variable
#'
#'\eqn{L_{n}(\gamma_{l};\lambda_{n})=||X_{l}-U\lambda_{l}||_{2}^{2}+\lambda_{n}\sum_{j=1}^{p} \omega_{njl} ||\gamma_{jl}||_{2}},where \eqn{\omega_{jnl}=||\gamma_{jl}||_{2}^{-1}},if \eqn{||\gamma_{jl}||_{2}>0},\eqn{\omega_{jnl}=infty},if \eqn{||\gamma_{jl}||_{2}=0}
#'
#'By minimizing the penalized objective function with a group Lasso penalty we by minimizing the penalized objective function with a group Lasso penalty. And then we use the selected IV for \eqn{\beta} in the model with two-stage least squares (2SLS).
#'
#' @return An object of type \code{naivereg} which is a list with the following
#' components:
#' \item{coefficients}{a vector of coefficients}
#' \item{ste}{the standard deviation of the coefficients}
#' \item{n}{number of samples}
#' \item{degree}{degree of B-splines}
#' \item{criterion}{The criterion by which to select the regularization parameter. One of "AIC", "BIC", "GCV", "AICc", or "EBIC"; default is "BIC".}
#' \item{ind}{the index of selected instrument variables}
#' \item{ind.b}{the index of selected instrument variables after B-splines}
#' \item{res}{the difference between the predicted y and the actual y}
#' @author Qingliang Fan, KongYu He, Wei Zhong
#' @references Q. Fan and W. Zhong (2017), “Nonparametric Additive Instrumental Variable Estimator: A Group Shrinkage Estimation Perspective,” Journal of Business & Economic Statistics, doi: 10.1080/07350015.2016.1180991.
#' @references Caner, M. and Fan, Q. (2015), Hybrid GEL Estimators: Instrument Selection with Adaptive Lasso, Journal of Econometrics, Volume 187, 256–274.
#' @examples
#'#naive regression
#'library(naivereg)
#'data("naivedata")
#'x=naivedata[,1]
#'y=naivedata[,2]
#'z=naivedata[,3:102]
#'#estimate with intercept
#'naive_intercept= naivereg(y,x,z)
#'#estimate without intercept,criterion:EBIC
#'naive_without_intercept = naivereg(y,x,z,intercept=FALSE,criterion='EBIC')
#' @importFrom splines bs
#' @importFrom grpreg grpreg select
#' @importFrom gmm gmm gel
#' @importFrom stats AIC BIC logLik
#' @export


naivereg <- function(y,x,z,max.degree=10,intercept=TRUE,criterion=c("BIC","AIC","GCV","AICc","EBIC"),
                       df.method=c("default","active"),
                       penalty = c("grLasso", "grMCP", "grSCAD", "gel", "cMCP"),endogenous.index=c(),IV.intercept=FALSE,family=c("gaussian", "binomial", "poisson")) {
  if(dim(as.matrix(z))[2]<dim(as.matrix(x))[2]){
    stop('The number of selected instruments is less than the number of endogenous variables, please find other instruments')
  }
  x=as.matrix(x)
  z=as.matrix(z)
  y=as.matrix(y)
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
    endogenous.index=rep(1,dim(as.matrix(x))[2])
  }
  if(missing(family)){
    family='gaussian'
  }

  z=cbind(z,x[,endogenous.index!=1])
  x=x[,endogenous.index]


  IV<-IVselect(z,x,criterion=criterion,df.method=df.method,penalty=penalty,IV.intercept=IV.intercept,family=family)
  naivelasso <- est.iv(y,x,IV$IVselect,intercept)

  list(coefficients=naivelasso$coefficients,ste=naivelasso$ste,n=naivelasso$n,degree=naivelasso$degree,criterion=criterion,
       ind=IV$ind,ind.b=IV$ind.b,res=naivelasso$residuals)

}


