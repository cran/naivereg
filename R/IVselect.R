#' @title Selecting instrument variables using group lasso and B-splines
#' @description Using group lasso and B-splines to obtain the group Lasso estimator where BIC or EBIC are applied to choose the tuning parameters  and (degree of B-splines) . And using the group
#' Lasso estimator to obtain the valid instrument variables
#' @param z The instrument variables matrix
#' @param x The design matrix, without an intercept
#' @param max.degree The upper limit value of degree of B-splines when using BIC/AIC to choose the tuning parameters, default is BIC.
#' @param criterion The criterion by which to select the regularization parameter. One of "AIC", "BIC", "GCV", "AICc", or "EBIC"; default is "BIC".
#' @param df.method How should effective model parameters be calculated? One of: "active", which counts the number of nonzero coefficients; or "default", which uses the calculated df returned by grpreg. default is "default".
#' @param penalty The penalty to be applied to the model. For group selection, one of grLasso, grMCP, or grSCAD. For bi-level selection, one of gel or cMCP. Default is " grLasso".
#' @param endogenous.index Specify which variables in design matrix are endogenous variables, the  variable corresponds to the value 1 is endogenous variables,  the  variable corresponds to the value 0  is exogenous variable, the default is all endogenous variables
#' @param IV.intercept Intercept of instrument variables, default is “FALSE”
#' @param family Either "gaussian" or "binomial", depending on the response.default is " gaussian "
#' @details See naivereg
#' @return An object of type \code{IVselect} which is a list with the following
#' components:
#' \item{degree}{degree of B-splines}
#' \item{criterion}{The criterion by which to select the regularization parameter. One of "AIC", "BIC", "GCV", "AICc", or "EBIC"; default is "BIC".}
#' \item{ind}{the index of selected instrument variables}
#' \item{ind.b}{the index of selected instrument variables after B-splines}
#' \item{IVselect}{The instrument variables after B-splines}
#' @author Qingliang Fan, KongYu He, Wei Zhong
#' @references Q. Fan and W. Zhong (2017), “Nonparametric Additive Instrumental Variable Estimator: A Group Shrinkage Estimation Perspective,” Journal of Business & Economic Statistics, doi: 10.1080/07350015.2016.1180991.
#' @references Caner, M. and Fan, Q. (2015), Hybrid GEL Estimators: Instrument Selection with Adaptive Lasso, Journal of Econometrics, Volume 187, 256–274.
#' @examples
#'#IV selecting with group Lasso an B-splines
#'library(naivereg)
#'data("naivedata")
#'x=naivedata[,1]
#'y=naivedata[,2]
#'z=naivedata[,3:102]
#'IV = IVselect(z,x)
#'IV$IVselect	#show the IV selected after B-splines
#' @export

IVselect<-function(z,x,max.degree=10,criterion=c("BIC","AIC","GCV","AICc","EBIC"),
                   df.method=c("default","active"),penalty = c("grLasso", "grMCP", "grSCAD", "gel", "cMCP"),
                   endogenous.index=c(),IV.intercept=FALSE,family=c("gaussian", "binomial", "poisson")){
  x=as.matrix(x)
  z=as.matrix(z)
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
  if(IV.intercept==TRUE){
    z <- cbind(matrix(1,nrow=dim(z)[1],ncol=1),z)
  }

  z=cbind(z,x[,endogenous.index!=1])
  x=x[,endogenous.index]

  degreeIC <- c(Inf,Inf)
  x <- as.matrix(x)
  z <- as.matrix(z)
  for(i in seq(3,max.degree,1)){
    IVbs=IVbs(z,i)
    fit=grpreg(IVbs$Z,x,IVbs$group,penalty=penalty)
    ll <- logLik(fit, df.method=df.method)
    df <- as.numeric(attr(ll,"df"))
    d <- dim(fit$beta)
    q <- if (length(d)==2) d[1]-1 else d[2]-1
    j <- if(fit$family=="gaussian") df-2 else df-1

    IC <- switch(criterion,
                 AIC = AIC(ll),
                 BIC = BIC(ll),
                 GCV = (1/fit$n) * (-2) * as.numeric(ll) / (1-df/fit$n)^2,
                 AICc = AIC(ll) + 2*df*(df+1)/(fit$n-df-1),
                 EBIC = BIC(ll) + 2*(lgamma(q+1) - lgamma(j+1) - lgamma(q-j+1)))
    degreeIC[i]<-min(IC)
  }
  degree <- which.min(degreeIC)
  IVbs <- IVbs(z,degree)
  p<-dim(as.matrix(z))[2]
  if(dim(x)[2]==1){
    fit<-grpreg(IVbs$Z,x,IVbs$group,penalty=penalty,family=family)

  ind.b <- (1:(p*degree))[(abs(select(fit,criterion)$beta[-1])>1e-5)]   # index of the selected B-splined variables with nonzero coefficients
  ind <- ind.b[seq(degree,p*degree,degree)]/degree
  ind <- ind[c(1:(length(ind.b)/degree))]
  }else
  {
    fit<-grpreg(IVbs$Z,x,IVbs$group,penalty=penalty,family=family)
    beta=matrix(0,dim(fit$beta[1,,])[1],dim(fit$beta[1,,])[2])
    for (h in seq(1,dim(x)[2])){
      beta=beta+fit$beta[h,,]^2
    }
    beta=sqrt(beta)
    ll <- logLik(fit, df.method=df.method)
    df <- as.numeric(attr(ll,"df"))
    d <- dim(fit$beta)
    q <- if (length(d)==2) d[1]-1 else d[2]-1
    j <- if(fit$family=="gaussian") df-2 else df-1

    IC <- switch(criterion,
                 AIC = AIC(ll),
                 BIC = BIC(ll),
                 GCV = (1/fit$n) * (-2) * as.numeric(ll) / (1-df/fit$n)^2,
                 AICc = AIC(ll) + 2*df*(df+1)/(fit$n-df-1),
                 EBIC = BIC(ll) + 2*(lgamma(q+1) - lgamma(j+1) - lgamma(q-j+1)))
    i <- which.min(IC)
    ind.b <- (1:(p*degree))[(abs(beta[,i][-1])>1e-5)]   # index of the selected B-splined variables with nonzero coefficients
    ind <- ind.b[seq(degree,p*degree,degree)]/degree
    ind <- ind[c(1:(length(ind.b)/degree))]
  }
    list(degree=degree,criterion=criterion,ind=ind,ind.b=ind.b,IVselect=(IVbs$Z)[,ind.b])
}
