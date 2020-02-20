\name{DSLIVEdata}
\alias{DSLIVEdata}
\title{The data generating for the DSLIVE}
\description{
##The data generation process is as follows
mu<-rep(0,200)

var<-matrix(,200,200)

for(i in 1:200){

  for(j in 1:200){

    var[i,j] <-0.5^(abs(i-j))

  }

}

x<-mvrnorm(100,mu,var)

mu<-rep(0,20)

varz<-matrix(,20,20)

for(i in 1:20){

  for(j in 1:20){

    varz[i,j] <-0.5^(abs(i-j))

  }

}

z<-mvrnorm(100,mu,varz)

mu1<-c(0,0)

v<-c(1,0.9,0.9,1)

var1<-matrix(v,2,2)

epsilon<-mvrnorm(100,mu1,var1)

D=rep(0,nrow(x))

p=.9*x[,4]+.8*x[,1]+1.96*x[,2]+1.85*x[,3]+.7*x[,5]+1.16*z[,1]+.95*z[,3]+1.7*z[,2]+epsilon[,2]

for(kk in 1:length(p)){

  D[kk]<-rbinom(n=1,size =  1,prob =  exp(p[kk])/(1+exp(p[kk])))

}

y=0.75*D+3*x[,1]+2*x[,5]+1.5*x[,4]+0*3*x[,7]+0*1.5*x[,8]+.15*x[,2]+.18*x[,3]+epsilon[,1]
  \itemize{
    \item Columns 1: Response variable y, an Nx1 vector.
    \item Columns 2-201: control variables x, an Nxp1 matrix.
    \item Columns 202-221: Instrumental variables,an Nxp2 matrix.
    \item Columns 222: Endogenous treatment variable,the value of endogenous variable is 0 or 1 (binary).
    }
}
\usage{data(DSLIVEdata)}
\keyword{datasets}

