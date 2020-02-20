#' @title DS-LIVE
#' @description Double selection plus logistic regression instrumental variable estimator (DS-LIVE). A three-step approach to estimate the dummy endogenous treatment effect using high-dimensional instruments in a penalized logistic regression model and double selection. This method accommodates the binary endogenous variable as well as the high-dimensionality for both the reduced form and structural equation models.
#' @param y Response variable, an N x 1 vector.
#' @param x Control variables, an N x p1 matrix.
#' @param z Instrumental variables, an N x p2 matrix.
#' @param D Endogenous treatment variable, the value of endogenous variable is 0 or 1 (binary).
#' @param criterion The criterion by which to select the regularization parameter. One of "BIC", "CV", CV means cross-validation, default is "BIC".
#' @param penalty This parameter takes effect when the creterion is CV. Quantitative for family="gaussian", or family="poisson" (non-negative counts). For family="binomial" should be either a factor with two levels, or a two-column matrix of counts or proportions (the second column is treated as the target class; for a factor, the last level in alphabetical order is the target class). For family="multinomial", can be a nc>=2 level factor, or a matrix with nc columns of counts or proportions. For either "binomial" or "multinomial", if y is presented as a vector, it will be coerced into a factor. For family="cox", y should be a two-column matrix with columns named 'time' and 'status'. The latter is a binary variable, with '1' indicating death, and '0' indicating right censored. The function Surv() in package survival produces such a matrix. For family="mgaussian", y is a matrix of quantitative responses.
#' @param family Only applied to the first step in the algorithm, the regression of y on x. Quantitative for family="gaussian", or family="poisson" (non-negative counts). For family="binomial" should be either a factor with two levels, or a two-column matrix of counts or proportions (the second column is treated as the target class; for a factor, the last level in alphabetical order is the target class). For family="multinomial", can be a nc>=2 level factor, or a matrix with nc columns of counts or proportions. For either "binomial" or "multinomial", if y is presented as a vector, it will be coerced into a factor. For family="cox", y should be a two-column matrix with columns named 'time' and 'status'. The latter is a binary variable, with '1' indicating death, and '0' indicating right censored. The function Surv() in package survival produces such a matrix. For family="mgaussian", y is a matrix of quantitative responses.
#' @param alpha Tuning parameter for the Mnet estimator which controls the relative contributions from the MCP/SCAD penalty and the ridge, or L2 penalty. alpha=1 is equivalent to MCP/SCAD penalty, while alpha=0 would be equivalent to ridge regression. However, alpha=0 is not supported; alpha may be arbitrarily small, but not exactly 0.
#' @param gamma The tuning parameter of the MCP/SCAD penalty. Default is 3.7.
#' @param nfolds This parameter takes effect when the creterion is CV. The response number of folds - default is 10. Although nfolds can be as large as the sample size (leave-one-out CV), it is not recommended for large datasets. Smallest value allowable is nfolds=3.
#' @param nlambda The number of lambda values, default is 100.
#' @param ... other arguments, see help(glmnet) or help(cv.ncvreg).
#' @details The DS-IV algorithm consists of the following three steps: In the first step, it estimates the coefficients (betaX) and select the important control variables set (denoted by c1) which are helpful to predict the outcome variable y using regularization methods for the data (y; x). In the second step, using a penalized logistic regression model, it selects both important control variables x (the selected control variables set is denoted by cx) and instrumental variables z for the endogenous treatment D. This step is crucial in the algorithm. Because it can estimate the optimal instrument using high-dimensional IVs as well as select additional important control variables which might be missed in the first step but are nonetheless important to the treatment variable. In the third step, it computes the post-double-selection LIVE estimator for the dummy endogenous treatment effect based on the predicted treatment variable D and the union of selected control variables in the first two variable selection steps denoted by c3 = (c1 union cx).
#' @return An object of type \code{DSLIVE} which is a list with the following
#' components:
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
#' @references Wei Zhong, Wei Zhou, Qingliang Fan and Yang Gao (2020), “Dummy Endogenous Treatment Effect Estimation Using High-Dimensional Instrumental Variables”, working paper.
#' @examples
#' library(naivereg)
#' data("DSLIVEdata")
#' y=DSLIVEdata[,1]
#' x=DSLIVEdata[,2:201]
#' z=DSLIVEdata[,202:221]
#' D=DSLIVEdata[,222]
#' res = DSLIVE(y,x,z,D,family='gaussian', criterion='BIC')
#' res$c1 # Variable indication of the selected in the first step (control variables x).
#' res$cx # Variable indication of selected control variables in the second step.
#' res$cz # Variable indication of selected instrumental variables in the second step.
#' res$c3 # Union of c1 and cx on control variables.
#' @export
######################################
DSLIVE <- function(y,x,z,D,criterion=c("BIC","CV"),penalty=c("SCAD","MCP","lasso"),family = c("gaussian", "binomial", "poisson", "multinomial","cox", "mgaussian"),alpha = 1
                 , gamma = 3.7,nfolds=10,nlambda = 100,...) {
  if(missing(criterion)){
    criterion='BIC'
  }
  if(missing(penalty)){
    penalty='SCAD'
  }
  if(missing(family)){
    family='gaussian'
  }

  #the first step
  data<-as.data.frame(cbind(y,x))
  x1=model.matrix(y~.,data)[,-1]
  y1=data[,1]
  if(criterion=='CV'){
    cv.out=cv.ncvreg(x1,y1,nfolds=nfolds,family=family,penalty=penalty,...)
    best.lambda<-cv.out$lambda.min
    coef2<-ncvreg(x1,y1,family = family,penalty = penalty,gamma = gamma,alpha = alpha,lambda=best.lambda,...)
    coef.scad<-coef2$beta[-1]
    c1<-which(as.numeric(coef.scad)!=0)
  }else{
    coef2<-glmnet(x1,y1,family=family,alpha=alpha,...)
    ###BIC
    pr<-predict(coef2,x1, s = NULL,type="link", exact = FALSE)
    loss<-colMeans((pr-y1)^2)
    BIC1<-rep(0,ncol(coef2$beta))
    for (n in 1:ncol(coef2$beta)) {
      BIC1[n]<-log(loss[n])+coef2$df[n]*(log(n)+0.1*log(ncol(x1)))/n#
    }
    #####
    best.lambda<-coef2$lambda[which.min(BIC1)]#
    coef2<-glmnet(x1,y1,family=family,alpha=alpha,lambda=best.lambda,...)
    c1<-which(as.numeric(coef2$beta)!=0)
  }

  #the second step
  data<-as.data.frame(cbind(D,x,z))
  x2=model.matrix(D~.,data)[,-1]
  y2=data[,1]
  if(criterion=='CV'){
    cv.out=cv.ncvreg(x2,y2,nfolds=nfolds,family="binomial",penalty=penalty,...)
    best.lambda<-cv.out$lambda.min
    coef2<-ncvreg(x2,y2,family = "binomial",penalty = penalty,gamma = gamma,alpha = alpha,lambda=best.lambda,...)
    coef.scad<-coef2$beta[-1]
    c2<-which(as.numeric(coef.scad)!=0)
  }else{
    ffit<-glmnet(x2,y2,family="binomial",alpha=alpha,nlambda=nlambda,lambda.min.ratio=0.01,...)
    ###  compute  BIC
    pr1<-predict(ffit,x2, s = NULL,type="response", exact = FALSE)
    for (h in 1:nrow(pr1)) {
      if(D[h]==0){
        pr1[h,]=1-pr1[h,]
      }
    }
    loss1<-colSums(log(pr1))
    BIC2<-rep(0,ncol(ffit$beta))
    for (nn in 1:ncol(ffit$beta)) {
      BIC2[nn]<--2*loss1[nn]+ffit$df[nn]*log(nrow(x2))
    }
    #####
    best.lambda1<-ffit$lambda[which.min(BIC2)]# optimal  BIC
    ffit<-glmnet(x2,y2,family="binomial",alpha=alpha,lambda=best.lambda1,...)
    c2<-which(as.numeric(ffit$beta)!=0)
  }


  #refit hat D
  se.va.scad<-as.matrix(x2[,c2])#which(coef.scad!=0)
  num.s<-ncol(se.va.scad)
  glm.scad=glm(D~se.va.scad[,1:num.s],family=binomial("logit"),control=list(maxit=500))
  DDD=glm.scad$fitted.value
  x3<-cbind(x,z)
  cx<-c2[c2<=ncol(x1)]
  cx<-cx[cx!=0]
  cz<-c2[c2>ncol(x1)]

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
  ctemp<-c(c1,cx)
  c3<-sort(ctemp[ctemp>0])  #
  #####
  xxx<-x[,c3]#final control variables

  datanew1<-as.data.frame(cbind(y,DDD,xxx))
  fitnew1<-lm(y~.-1,datanew1)
  beta<-fitnew1$coef
  betaD<-beta[1]
  betaX<-beta[-1]
  list(betaD=betaD,betaX=betaX,c1=c1,cx=cx,cz=cz,c2=c2,c3=c3,family = family,criterion=criterion)
}
