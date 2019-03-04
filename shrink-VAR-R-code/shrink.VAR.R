# 2007 Rainer Opgen-Rhein
# License: GNU GPL 2 or later


shrink.VAR<-function(data, lambda, lambda.var){

# split data in past and future
  p<-ncol(data)

  X<-data
  if(is.longitudinal(X)==FALSE){
    stop("input is not a longitudinal object")
  }
  n.v<-ncol(X)
  rep<-get.time.repeats(X)$repeats[1]
  Xp<-(X[1:(nrow(X)-rep),])
  Xf<-(X[-(1:rep),])
  XpXf<-cbind(Xp,Xf)


# calculate beta-matrix

covXpXf<-cov.shrink(XpXf, lambda=lambda, lambda.var=lambda.var)


l.est<-attr(covXpXf, "lambda")
l.var.est<-attr(covXpXf, "lambda.var")

B_est<-solve(covXpXf[1:p,1:p])%*%covXpXf[1:p,(p+1):(2*p)]

print("calculate partial correlations")
##calculate partial correlations

## use B_star
B.star<-matrix(ncol=p,nrow=p)

 for(i in 1:p){
  b.star<-mvr.shrink(Xf[,i],Xp,lambda=l.est, lambda.var=l.var.est,verbose=FALSE)[,-1]
  B.star[,i]<-b.star
  print(i)
 } 
  r.mat.abs<-sqrt(B_est*B.star)
  r.mat<-sign(B_est)*r.mat.abs
 
 
results<-list(B_est=B_est, r.mat=r.mat, lambda=l.est, lambda.var=l.var.est)
return(results)
}


test.VAR<-function(out.Var, df=10, plot.locfdr=1){

B_est<-out.Var$B_est
r.mat<-out.Var$r.mat
p<-ncol(B_est)


#### list nodes

    value<-as.vector(r.mat)
    leng<-nrow(r.mat)
    indexes <- matrix(nrow=leng^2,ncol=2)
    indexes[,2]<-rep(1:leng,each=leng)
    indexes[,1]<-rep(1:leng,leng)
    colnames(indexes) <- c("node1_past", "node2_future")
    node.list <- cbind(value, indexes)
    sort.idx <- order(-abs(node.list[, 1]))
    node.list <- node.list[sort.idx, ]
    

## calculate significance

    pc <- node.list[, 1]

# locfdr
     require("fdrtool")
     if (any(pc > 1) || any(pc < -1))
        stop("Data out of range: input correlations must be in [-1; 1]")
   out <-fdrtool(z.transform(pc),"correlation")   # 

   prob.nonzero <- 1 - out$lfdr
   node.sig <- cbind(node.list, out$pval, out$qval, prob.nonzero)#

    node.sig<-as.data.frame(node.sig)

return(node.sig)
}   


