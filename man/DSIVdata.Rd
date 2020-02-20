\name{DSIVdata}
\alias{DSIVdata}
\title{The data generating for the DSIV}
\description{
##The data generation process is as follows

library(MASS)

mu<-rep(0,50)

var<-matrix(,50,50)

for(i in 1:50){

  for(j in 1:50){

    var[i,j] <-0.5^(abs(i-j))

  }

}

x<-mvrnorm(n,mu,var)#generate x

mu<-rep(0,20)

varz<-matrix(,20,20)

for(i in 1:20){

  for(j in 1:20){

    varz[i,j] <-0.5^(abs(i-j))

  }

}

z<-mvrnorm(n,mu,varz)#generate iv

mu1<-c(0,0)

v<-c(1,0.9,0.9,1)

var1<-matrix(v,2,2)

epsilon<-mvrnorm(n,mu1,var1)#generate error term

D=1.9*x[,2]+2.5*x[,3]+1.4*x[,5]+x[,6]+x[,1]+1.6*z[,1]+1.9*z[,3]+1.7*z[,2]+epsilon[
,2]

y=0.75*D+1*x[,1]+2*x[,6]+.11*x[,2]+.18*x[,3]+.12*x[,5]+epsilon[,1]
  \itemize{
    \item Columns 1: Response variable y, an Nx1 vector.
    \item Columns 2-51: control variables x, an Nxp1 matrix.
    \item Columns 52-71: Instrumental variables,an Nxp2 matrix.
    \item Columns 72: Endogenous treatment variable.
    }
}
\usage{data(DSIVdata)}
\keyword{datasets}

