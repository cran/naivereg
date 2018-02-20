
est.iv <- function (y, x, z,intercept) {
  x <- as.matrix(x)
  n <- length(y)
  if(intercept==TRUE){
    x <- cbind(matrix(1,nrow=n,ncol=1),x)
  }
  p <- ncol(x)
  z <- as.matrix(z)
  invZtZ <- solve(crossprod(z))
  XtZ <- crossprod(x, z)
  V <- chol2inv(chol(XtZ %*% invZtZ %*% t(XtZ)+diag(1e-6,p)))
  b <- V %*% XtZ %*% invZtZ %*% crossprod(z, y)
  residuals <- y - x %*% b
  s2 <- sum(residuals^2)/(n - p)
  V <- s2*V
  ste <- sqrt(diag(V))
  list(coefficients=b,ste=ste,n=n,residuals=as.vector(residuals),instruments=z)
}
