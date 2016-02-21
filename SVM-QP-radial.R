rm(list=ls(all=TRUE))
library(e1071)
library(quadprog)
library(lpSolve)

std<-function(x) { return((x-mean(x))/sd(x)) }
soft<-function(x, lam=0) { return(ifelse(x>lam,x-lam,0)) }
sign<-function(x) { return(ifelse(x>=0,1,-1)) }
K<-function(x,y,kernel="radial",sigma2=ifelse(is.vector(x),1,ncol(x))) {
  if(is.vector(x)) {
    if (kernel=="linear") { return ( sum(x*y) ) } else { return ( exp(-sum((x-y)^2)/sigma2) ) }
  } else {
    if (kernel=="linear") return ( x %*% t(y) ) else {
      res<-matrix(0,nrow(x),nrow(y))
      for ( i in 1:nrow(x) ) {
        for ( j in 1:nrow(y) ) { 
          res[i,j]<-exp(-sum((x[i,]-y[j,])^2)/sigma2) 
        }
      }
      return (res)
    }
  }
}


### kernel SVM ###
### penalize intercept ###

inf.small=1e-3
svm.prog.ker=function(xx,yy,lambda=1,kernel="radial"){
  n=nrow(xx)
  p=ncol(xx)
  
  # ker.mat=outer(1:n,1:n,FUN=Vectorize(function(i,j) K(xx[i,],xx[j,],kernel)))
  ker.mat=K(xx,xx,kernel)
  Dmat= (ker.mat+matrix(1,n,n))*(yy%*%t(yy))+diag(inf.small,n)
  dvec=rep(1,n)
  Amat=cbind(diag(1,n),diag(-1,n))
  bvec=c(rep(0,n),rep(-1/lambda,n))
  sol=solve.QP(Dmat,dvec,Amat,bvec,meq=0)
  alpha=sol$solution
  normsq=sum(diag(alpha)%*%Dmat%*%diag(alpha))
  
  list(alpha=alpha,normsq=normsq)
}

########################The simulation#################
seed_num = 12

set.seed(seed_num)
nrow = 500 #put number of rows here to change
dataset = 3
tuneP = 2
per = 1

###########################################################
if (dataset==1){
  x=matrix(runif(nrow*10,-1,1),nrow)
  y=sign(x[,1])
  y.ind=sample(1:nrow,nrow/200,replace=F)
  y[y.ind]=-y[y.ind]
  
  x.test=x[1:(nrow/2),]
  y.test=y[1:(nrow/2)]
  x=x[-(1:(nrow/2)),]
  y=y[-(1:(nrow/2))]
}

###########################################################
if (dataset==2){
  n=nrow/4
  sigma <- 2.5
  x1 <- rnorm(n,-4,sigma)
  y1 <- rnorm(n,-2,sigma)
  x2 <- rnorm(n,2,sigma)
  y2 <-  rnorm(n,4,sigma)
  x=c(x1,x2,y1,y2)
  x=matrix(x,ncol=2)
  y=c(-rep(1,n),rep(1,n))
  
  x1t <- rnorm(n,-4,sigma)
  y1t <- rnorm(n,-2,sigma)
  x2t <- rnorm(n,2,sigma)
  y2t <-  rnorm(n,4,sigma)
  tempx=c(x1t,x2t,y1t,y2t)
  x.test=matrix(tempx,ncol=2)
  y.test=c(-rep(1,n),rep(1,n))
}
###########################################################
if (dataset==3){
  n=nrow/4
  r1=3;
  r2=5;
  r3=7;
  angle1 = runif(n,0,2*pi);
  angle2 = runif(n,0,2*pi);
  x1 <- r1+(r2-r1)*cos(angle1)*rnorm(n,1,0.2)
  y1 <- r1+(r2-r1)*sin(angle1)*rnorm(n,1,0.2)
  x2 <- r1+(r3-r1)*cos(angle2)*rnorm(n,1,0.2)
  y2 <- r1+(r3-r1)*sin(angle2)*rnorm(n,1,0.2)
  x=c(x1,x2,y1,y2)
  x=matrix(x,ncol=2)
  y=c(-rep(1,n),rep(1,n))
  
  angle3 = runif(n,0,2*pi);
  angle4 = runif(n,0,2*pi);
  x1t <- r1+(r2-r1)*cos(angle3)*rnorm(n,1,0.2)
  y1t <- r1+(r2-r1)*sin(angle3)*rnorm(n,1,0.2)
  x2t <- r1+(r3-r1)*cos(angle4)*rnorm(n,1,0.2)
  y2t <- r1+(r3-r1)*sin(angle4)*rnorm(n,1,0.2)
  tempx=c(x1t,x2t,y1t,y2t)
  x.test=matrix(tempx,ncol=2)
  y.test=c(-rep(1,n),rep(1,n))
}






########################The simulation#################
times=15
seed_num = sample(1:100,times)
type = "radial"
k1 = 10
k2 = 2

error0=0
error1=0
error2=0
error3=0
iter = 25
t_total = c(0,0,0,0)
t_indiv = c(0,0,0,0)

ker.mat=matrix(0,nrow,nrow)

for (tim in 1:times) {
  ############# |yf(x)|+random ############
  set.seed(seed_num[tim])
  samp.set=sample(1:(nrow/2),20,F)
  cand.set=setdiff(1:(nrow/2),samp.set)
  
  test.err1=NULL
  ptm11=0
  ptm12=0
  ptm13=0
  for(i in 1:iter) {
    ptm_temp = proc.time()[1]
    xx=x[samp.set,]
    yy=y[samp.set]
    m1=svm.prog.ker(xx,yy,lambda=1)  
    alpha=m1$alpha
    c.normsq=m1$normsq

    if(i==1) {
      sel.set2=samp.set
      candQ.set = cand.set
    }
    ker.mat[samp.set,candQ.set]=K(x[samp.set,],x[candQ.set,],type)
    ptm11 <- ptm11 + proc.time()[1]-ptm_temp
    ptm_temp1 = proc.time()[1]
    
    ker.mat[sel.set2,(nrow/2+1):nrow]=K(x[sel.set2,],x.test,type)
    m.pred=sign(sum(alpha*yy)+apply(alpha*yy*ker.mat[samp.set,(nrow/2+1):nrow],2,sum))
    test.err1=c(test.err1,sum(y.test!=m.pred)/length(y.test))
    ptm12 <- ptm12 + proc.time()[1]-ptm_temp1
    ptm_temp2 = proc.time()[1]
    
    cd = length(cand.set)
    candQ.set=sample(cand.set,cd*per,F)
    cand.fit=sum(alpha*yy)+apply(alpha*yy*ker.mat[samp.set,candQ.set],2,sum)
    sel.set1=candQ.set[order(abs(y[candQ.set]*cand.fit))[1:k1]]
    sel.set2=sample(sel.set1,k2)
    samp.set=c(samp.set,sel.set2)
    cand.set=setdiff(cand.set,sel.set2)
    ptm13 <- ptm13 + proc.time()[1]-ptm_temp2
  }
  
  ############# |yf(x)| ############
  set.seed(seed_num[tim])
  samp.set=sample(1:(nrow/2),20,F)
  cand.set=setdiff(1:(nrow/2),samp.set)
  
  test.err2=NULL
  ptm21=0
  ptm22=0
  ptm23=0
  
  for(i in 1:iter) {
    ptm_temp = proc.time()[1]
    xx=x[samp.set,]
    yy=y[samp.set]
    m1=svm.prog.ker(xx,yy,lambda=1)  
    alpha=m1$alpha
    c.normsq=m1$normsq
    
    if(i==1){
      sel.set=samp.set 
      candQ.set = cand.set
    }
    ker.mat[samp.set,candQ.set]=K(x[samp.set,],x[candQ.set,],type)
    ptm21 <- ptm21 + proc.time()[1]-ptm_temp
    ptm_temp1 = proc.time()[1]
    
    ker.mat[sel.set,(nrow/2+1):nrow]=K(x[sel.set,],x.test,type)
    m.pred=sign(sum(alpha*yy)+apply(alpha*yy*ker.mat[samp.set,(nrow/2+1):nrow],2,sum))    
    test.err2=c(test.err2,sum(y.test!=m.pred)/length(y.test))
    ptm22 <- ptm22 + proc.time()[1]-ptm_temp1
    ptm_temp2 = proc.time()[1] 
    
    cd = length(cand.set)
    candQ.set=sample(cand.set,cd*per,F)
    cand.fit=sum(alpha*yy)+apply(alpha*yy*ker.mat[samp.set,candQ.set],2,sum)
    sel.set=candQ.set[order(abs(y[candQ.set]*cand.fit))[1:k2]]    
    samp.set=c(samp.set,sel.set)
    cand.set=setdiff(cand.set,sel.set)
    ptm23 <- ptm23 + proc.time()[1]-ptm_temp2
  }
  
  ############# random ############
  set.seed(seed_num[tim])
  samp.set=sample(1:(nrow/2),20,F)
  cand.set=setdiff(1:(nrow/2),samp.set)
  
  test.err0=NULL
  ptm01=0
  ptm02=0
  ptm03=0
  
  for(i in 1:iter) {
    ptm_temp = proc.time()[1]
    xx=x[samp.set,]
    yy=y[samp.set]
    m1=svm.prog.ker(xx,yy,lambda=1)  
    alpha=m1$alpha
    c.normsq=m1$normsq
    
    if(i==1) sel.set=samp.set 
    ker.mat[samp.set,cand.set]=K(x[samp.set,],x[cand.set,],type)
    ptm01 <- proc.time()[1]-ptm_temp
    ptm_temp1 = proc.time()[1]
    
    ker.mat[sel.set,(nrow/2+1):nrow]=K(x[sel.set,],x.test,type)
    m.pred=sign(sum(alpha*yy)+apply(alpha*yy*ker.mat[samp.set,(nrow/2+1):nrow],2,sum))
    test.err0=c(test.err0,sum(y.test!=m.pred)/length(y.test))
    ptm02 <- proc.time()[1]-ptm_temp1
    ptm_temp2 = proc.time()[1]
    
    
    sel.set=sample(cand.set,k2,F)
    samp.set=c(samp.set,sel.set)
    cand.set=setdiff(cand.set,sel.set)
    ptm03 <- proc.time()[1]-ptm_temp2
  }
  
  error0 = error0+test.err0
  error1 = error1+test.err1
  error2 = error2+test.err2
  
  t_indiv = rbind(c(ptm11+ptm13,ptm21+ptm23,ptm01+ptm03,0),t_indiv)
  
  t_total = t_total + c(ptm11+ptm13,ptm21+ptm23,ptm01+ptm03,0)
}

ptm_all1 = proc.time()[1]
svm(x,as.factor(y),cost=1,kernel=type)
ptm_all2 = proc.time()[1] 

t_total = t_total + c(0, 0, 0, (ptm_all2-ptm_all1)*iter)
t_aver = t_total/iter

#########
error0 = error0/times
error1 = error1/times
error2 = error2/times

plot(error0,type="l",lty=2,ylim=c(0,0.2))
lines(error1,lty=1,col=2) #red yf(x)+random
lines(error2,lty=1,col=3) #green yf(x)


