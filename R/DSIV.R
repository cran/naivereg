#' @title Double-Selection Plus Instrumental Variable Estimator
#' @description A three-step approach to estimate the endogenous treatment effect using high-dimensional instruments and double selection. It is applicable in the following scenarios: first, there is a known endogeneity problem for the treatment variable. Second, the treatment effect model has a large number of control variables, such as the large micro survey data.
#' @param y Response variable, an N x 1 vector.
#' @param x Control variables, an N x p1 matrix.
#' @param z Instrumental variables, an N x p2 matrix.
#' @param D Endogenous treatment variable.
#' @param family Quantitative for family="gaussian", or family="poisson" (non-negative counts). For family="binomial" should be either a factor with two levels, or a two-column matrix of counts or proportions (the second column is treated as the target class; for a factor, the last level in alphabetical order is the target class). For family="multinomial", can be a nc>=2 level factor, or a matrix with nc columns of counts or proportions. For either "binomial" or "multinomial", if y is presented as a vector, it will be coerced into a factor. For family="cox", y should be a two-column matrix with columns named 'time' and 'status'. The latter is a binary variable, with '1' indicating death, and '0' indicating right censored. The function Surv() in package survival produces such a matrix. For family="mgaussian", y is a matrix of quantitative responses.
#' @param criterion The criterion by which to select the regularization parameter. One of "BIC", "EBIC", default is "BIC".
#' @param alpha The elasticnet mixing parameter, with 0<=alpha<= 1. alpha=1 is the lasso penalty, and alpha=0 the ridge penalty.
#' @param nlambda The number of lambda values, default is 100.
#' @param ... other arguments, see help(glmnet).
#' @details The DS-IV algorithm consists of the following three steps: In the first step,
#'regress the outcome variable y on control variables x using the
#'regularization method, estimate the coefficients beta and select the important control
#'variables set denoted by c1. In the second step, regress the treatment variable
#'d on instrumental variables w and control variables x, estimate the
#'optimal instrument d and obtain the second important control variables set
#'denoted by cx. In the third step, obtain the DS-IV estimator of the endogenous
#'of the endogenous treatment effect based on the estimated optimal instrument d
#'and the union (c3) of the selected control variables.
#' @return An object of type \code{DSIV} which is a list with the following
#' components:
#' \item{yhat}{The estimated value of y.}
#' \item{betaD}{The coefficient of endogenous variable D.}
#' \item{betaX}{The coefficient of control variables x.}
#' \item{c1}{Variable indication of the selected in the first step (control variables x).}
#' \item{cx}{Variable indication of selected control variables in the second step.}
#' \item{cz}{Variable indication of selected instrumental variables in the second step.}
#' \item{c2}{Variable indication of the selected in the second step. The number less than or equal to p1 is an indication of control variables, the number greater than p1 and less than or equal to (p1 + p2) is an indication of instrument variables.}
#' \item{c3}{Union of c1 and cx on control variables.}
#' \item{family}{Same as above.}
#' \item{criterion}{Same as above.}
#' @author Qingliang Fan, KongYu He, Wei Zhong
#' @references Wei Zhong, Yang Gao, Wei Zhou and Qingliang Fan (2020), “Endogenous Treatment Effect Estimation Using High-Dimensional Instruments and Double Selection”, working paper
#' @examples
#'library(naivereg)
#'data("DSIVdata")
#'y=DSIVdata[,1]
#'x=DSIVdata[,2:51]
#'z=DSIVdata[,52:71]
#'D=DSIVdata[,72]
#'res = DSIV(y,x,z,D,family='gaussian', criterion='EBIC')
#'res$c1 #Variable indication of the selected in the first step (control variables x).
#'res$cx #Variable indication of selected control variables in the second step.
#'res$cz #Variable indication of selected instrumental variables in the second step.
#'res$c3 #Union of c1 and cx on control variables
#' @export
######################################
DSIV <- function(y,x,z,D,family = c("gaussian", "binomial", "poisson", "multinomial","cox", "mgaussian"),criterion=c("BIC","EBIC"),alpha = 1
                ,nlambda = 100,...) {
  if(missing(criterion)){
    criterion='EBIC'
  }
  if(missing(family)){
    family='gaussian'
  }

  #The first stage of variable selection of double selection method
  n = length(y)
  data<-as.data.frame(cbind(y,x))
  x1=model.matrix(y~.,data)[,-1]
  y1=data[,1]

  fit1<-glmnet(x1,y1,family=family,alpha=alpha,nlambda=nlambda,...)
  pr1<-predict(fit1,x1, s = NULL,type="link", exact = FALSE)
  loss<-colMeans((pr1-y1)^2)
  #calculate IC
  IC<-rep(0,ncol(fit1$beta))
  for (nn in 1:ncol(fit1$beta)) {
    IC[nn] <- switch(criterion,
                 BIC = n*log(loss[nn])+fit1$df[nn]*log(n),
                 EBIC = log(loss[nn])+fit1$df[nn]*(log(n)+0.2*log(ncol(x1)))/n)
  }

  best.lambda<-fit1$lambda[which.min(IC)]
  fit11<-glmnet(x1,y1,family=family,alpha=alpha,lambda=best.lambda,nlambda=nlambda,...)
  c1<-which(as.numeric(fit11$beta)!=0)

  #the second stage of variable selection of double  selection method
  data<-as.data.frame(cbind(D,x,z))
  x2=model.matrix(D~.,data)[,-1]
  y2=data[,1]
  fit2<-glmnet(x2,y2,family=family,alpha=alpha,nlambda=nlambda,...)
  pr2<-predict(fit2,x2, s = NULL,type="link", exact = FALSE)
  loss1<-colMeans((pr2-y2)^2)
  #calculate IC
  IC1<-rep(0,ncol(fit2$beta))
  for (nn in 1:ncol(fit1$beta)) {
    IC1[nn] <- switch(criterion,
                     BIC = n*log(loss1[nn])+fit2$df[nn]*log(n),
                     EBIC = log(loss1[nn])+fit2$df[nn]*(log(n)+0.2*log(ncol(x2)))/n)
  }

  best.lambda1<-fit2$lambda[which.min(IC1)]
  fit21<-glmnet(x2,y2,family=family,alpha=alpha,lambda=best.lambda1,nlambda=nlambda,...)

  x3<-cbind(x,z)
  c2<-which(as.numeric(fit21$beta)!=0)

  xxxx<-x3[,c2]
  ########refit and estimate the optimal IV
  datanew2<-as.data.frame(cbind(D,xxxx))

  fitnew2<-lm(D~.-1,datanew2)
  DDD<-fitnew2$fitted.values

  cx<-c2[c2<=ncol(x)]
  cx<-cx[cx!=0]
  cz<-c2[c2>ncol(x)]
  c1<-c1[c1<ncol(x)+1]

  for (i in 1:length(c1)) {
    for (j in 1:length(cx)) {
      if(is.na(cx[j])==T){
        cx<-0
        break
      }else if(cx[j]==c1[i]){
        cx[j]<-0
      }
    }
  }
  c3<-c(c1,cx)
  c4<-c3[c3>0]
  xx<-x[,c4]

  ###DS-IV estimator
  pp<-xx%*%solve(t(xx)%*%xx)%*%t(xx)
  ii<-matrix(rep(0,nrow(xx)^2),nrow(xx),nrow(xx))
  for(t in 1:nrow(xx)){
    ii[t,t]<-1
  }
  m<-ii-pp
  betaD<-solve(t(DDD)%*%m%*%D)%*%t(DDD)%*%m%*%y
  datanew<-as.data.frame(cbind(y,DDD,xx))
  fitnew<-lm(y~.-1,datanew)
  betaX <- fitnew$coef[-1]
  yhat<-as.vector(betaD)*D-xx%*%betaX
  list(yhat=yhat,betaD=betaD,betaX=betaX,c1=c1,cx=cx,cz=cz,c2=c2,c3=c4,family = family,criterion=criterion)
}
