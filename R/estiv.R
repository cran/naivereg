
est.iv <- function (y,x_endogenous,x_exogenous,xhat,intercept) {

  xhat_all =cbind(xhat,x_exogenous)
  x_all = cbind(x_endogenous,x_exogenous)
  x_all <- as.matrix(x_all)
  n <- length(y)

  if(intercept==TRUE){
    x_all <- cbind(x_all,matrix(1,nrow=n,ncol=1))
    xhat_all <- cbind(xhat_all,matrix(1,nrow=n,ncol=1))
  }
  p <- ncol(x_all)

  V = chol2inv(chol(t(xhat_all)%*%x_all+diag(1e-6,p)))
  b <- V%*%t(xhat_all)%*%y

  residuals <- y - x_all %*% b
  s2 <- sum(residuals^2)/(n - p)
  V <- s2*V
  ste <- sqrt(diag(V))

  list(coefficients=b,ste=ste,n=n,residuals=as.vector(residuals))
}
