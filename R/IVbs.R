

IVbs<-function(Z,degree){ # Z is the matrix
  p <-dim(as.matrix(Z))[2]
  Z.bs<-NULL
  group<-c()
  for (j in 1:p) {
    Z.bs <-cbind(Z.bs,bs(Z[,j], degree))
    group<-c(group,rep(j,degree))
  }
  list(Z.bs=Z.bs,group=group)
}
