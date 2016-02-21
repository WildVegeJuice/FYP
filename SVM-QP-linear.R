rm(list=ls(all=TRUE))
library(e1071)
library(quadprog)
library(lpSolve)

############### SVM ######################
## NO intercept
## ONLY linear SVM

svm.prog=function(xx,yy,lambda=1){
  n=nrow(xx)
  p=ncol(xx)
  
  Dmat=((xx*yy) %*% t(xx*yy)) +diag(1e-4,n)
  dvec=rep(1,n)
  Amat=cbind(diag(1,n),diag(-1,n))
  bvec=c(rep(0,n),rep(-1/lambda,n))
  sol=solve.QP(Dmat,dvec,Amat,bvec,meq=0)
  alpha=sol$solution
  beta=apply(alpha*yy*xx,2,sum)
  
  return(beta)
}


########################The simulation#################
set.seed(11)
nrow = 1000 #put number of rows here to change, data size is nrow/2
dataset = 2
tuneP = 2
times = 10
type = "linear"
ap3 = "kancolle"
k1 = 10
k2 = 2
per = 1
iter = 20


###########################################################
if (dataset==1){
  x1=matrix(runif(nrow*5,-1,1),nrow/2)
  x1[,1]=-sign(x1[,1])*x1[,1]  
  x2=matrix(runif(nrow*5,-1,1),nrow/2)
  x2[,1]=sign(x2[,1])*x2[,1]
  x=rbind(x1,x2)
  y=c(-rep(1,nrow/2),rep(1,nrow/2))
  
  y.ind=sample(1:(nrow/2),nrow/400,replace=F)
  y[y.ind]=-y[y.ind]
  
  x.test=x[(nrow/4+1):(3*nrow/4),]
  y.test=y[(nrow/4+1):(3*nrow/4)]
  x=x[-((nrow/4+1):(3*nrow/4)),]
  y=y[-((nrow/4+1):(3*nrow/4))]
  ptm = proc.time()[1]
}

###########################################################
if (dataset==2){
  n=nrow/4
  sigma <- 2.2
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
  ptm = proc.time()[1]
}


########################The simulation#################
seed_num = sample(1:100,times)
error0=0; error1=0; error2=0; error3=0
e0=c(); e1=c(); e2=c(); e3=c()
t_total = c(0,0,0,0,0)
t_indiv = c(0,0,0,0,0)

for (tim in 1:times) {
  ############# |yf(x)|+random ############
  set.seed(seed_num[tim])
  samp.set=sample(1:(nrow/2),20,F)
  cand.set=setdiff(1:(nrow/2),samp.set)
  test.err1=NULL
  ptm11=0; ptm12=0; ptm13=0
  
  for(i in 1:iter) {
    ptm_temp = proc.time()[1]
    xx=cbind(1,x[samp.set,])
    yy=y[samp.set]
    beta=svm.prog(xx,yy,lambda=1)	
    ptm11 <- ptm11 + proc.time()[1]-ptm_temp
    ptm_temp1 = proc.time()[1]
    
    m.pred=sign(cbind(1,x.test)%*%beta)    
    test.err1=c(test.err1,sum(y.test!=m.pred)/length(y.test))
    ptm12 <- ptm12 + proc.time()[1]-ptm_temp1
    ptm_temp2 = proc.time()[1]
    
    cd = length(cand.set)
    candQ.set=sample(cand.set,cd*per,F)
    cand.fit=y[candQ.set]*(cbind(1,x[candQ.set,])%*%beta)

    rank.set=c(1:length(candQ.set))[cand.fit<1]
    if (length(rank.set) > k1){
      sel.set1=candQ.set[rank.set[order(abs(cand.fit[rank.set]))[1:k1]]]
    } else {
      sel.set1=candQ.set[order(abs(cand.fit))[1:k1]]
    }
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
  ptm21=0; ptm22=0; ptm23=0
  
  for(i in 1:iter) {
    ptm_temp = proc.time()[1]
    xx=cbind(1,x[samp.set,])
    yy=y[samp.set]
    beta=svm.prog(xx,yy,lambda=1)  
    ptm21 <- ptm21 + proc.time()[1]-ptm_temp
    ptm_temp1 = proc.time()[1]
    
    m.pred=sign(cbind(1,x.test)%*%beta)    
    test.err2=c(test.err2,sum(y.test!=m.pred)/length(y.test))
    ptm22 <- ptm22 + proc.time()[1]-ptm_temp1
    ptm_temp2 = proc.time()[1]  
    
    cd = length(cand.set)
    candQ.set=sample(cand.set,cd*per,F)
    cand.fit=y[candQ.set]*(cbind(1,x[candQ.set,])%*%beta)
    
    rank.set=c(1:length(candQ.set))[cand.fit<1]
    if (length(rank.set) > k2){
      sel.set=candQ.set[rank.set[order(abs(cand.fit[rank.set]))[1:k2]]]
    } else {
      sel.set=candQ.set[order(abs(cand.fit))[1:k2]]
    }   
    samp.set=c(samp.set,sel.set)
    cand.set=setdiff(cand.set,sel.set)
    ptm23 <- ptm23 + proc.time()[1]-ptm_temp2
  }
  
  ############# random ############
  set.seed(seed_num[tim])
  samp.set=sample(1:(nrow/2),20,F)
  cand.set=setdiff(1:(nrow/2),samp.set)
  ptm01=0; ptm02=0; ptm03=0
  test.err0=NULL
  
  for(i in 1:iter) {
    ptm_temp = proc.time()[1]
    xx=cbind(1,x[samp.set,])
    yy=y[samp.set]
    beta=svm.prog(xx,yy,lambda=1)       
    ptm01 <- proc.time()[1]-ptm_temp
    ptm_temp1 = proc.time()[1]
    
    m.pred=sign(cbind(1,x.test)%*%beta)
    test.err0=c(test.err0,sum(y.test!=m.pred)/length(y.test))
    ptm02 <- proc.time()[1]-ptm_temp1
    ptm_temp2 = proc.time()[1]
    
    sel.set=sample(cand.set,k2,F)
    samp.set=c(samp.set,sel.set)
    cand.set=setdiff(cand.set,sel.set)
    ptm03 <- proc.time()[1]-ptm_temp2
  }
  
  
  ############# |yf(x)|/K ############
  if (type=="lovelive"){
    set.seed(seed_num[tim])
    samp.set=sample(1:(nrow/2),20,F)
    cand.set=setdiff(1:(nrow/2),samp.set)
    test.err3=NULL
    ptm31=0; ptm32=0; ptm33=0
    
    diagK = cbind((rowSums(cbind(1,x)^2))^tuneP)
    
    for(i in 1:iter) {
      ptm_temp = proc.time()[1]
      xx=cbind(1,x[samp.set,])
      yy=y[samp.set]
      beta=svm.prog(xx,yy,lambda=1)  
      ptm_temp1 = proc.time()[1]
      ptm31 <- ptm31 + proc.time()[1]-ptm_temp  
      
      m.pred=sign(cbind(1,x.test)%*%beta)
      test.err3=c(test.err3,sum(y.test!=m.pred)/length(y.test))
      ptm_temp2 = proc.time()[1]  
      ptm32 <- ptm32 + proc.time()[1]-ptm_temp1
      
      cd = length(cand.set)
      candQ.set=sample(cand.set,cd*per,F)
      cand.fit=cbind(1,x[candQ.set,])%*%beta
      temp=diagK[candQ.set,]
      
      sel.set = candQ.set[order(abs(y[candQ.set]*cand.fit)/temp)[1:k2]]
      samp.set=c(samp.set,sel.set)
      cand.set=setdiff(cand.set,sel.set)
      ptm33 <- ptm33 + proc.time()[1]-ptm_temp2
    }
    
  }
  
  if (type=="lovelive"){
    error3 = error3+test.err3
    e3 = rbind(e3,test.err3)
  }
  
  error0 = error0+test.err0; error1 = error1+test.err1; error2 = error2+test.err2
  e0 = rbind(e0,test.err0,deparse.level=0); e1 = rbind(e1,test.err1,deparse.level=0); e2 = rbind(e2,test.err2,deparse.level=0)
  t_indiv = rbind(c(ptm11+ptm13,ptm21+ptm23,ptm01+ptm03,0,0),t_indiv)
  t_total = t_total + c(ptm11+ptm13,ptm21+ptm23, ptm01+ptm03,0,0)  
}

####################### direct solve ##################
ptm_all1 = proc.time()[1]
m0 = svm(x,as.factor(y),cost=1,kernel=type)
ptm_all2 = proc.time()[1] 

t_total = t_total + c(0, 0, 0, 0, (ptm_all2-ptm_all1)*iter)
t_aver = t_total/iter


######### Analyse error #########
se0=c(); se1=c(); se2=c(); se3=c()
for (i in 1:iter){
  se0 = c(se0,sd(e0[,i])/sqrt(times))
  se1 = c(se1,sd(e1[,i])/sqrt(times))
  se2 = c(se2,sd(e2[,i])/sqrt(times))
}

error0 = error0/times; error1 = error1/times; error2 = error2/times
if ( ap3 =="lovelive"){
  error3 = error3/times
  se3 = c(se3,sd(e3[,i]))
}

#########
#error0 = error0/times
#error1 = error1/times
#error2 = error2/times
#if (type=="lovelive"){
#  error3 = error3/times
#}
#plot(error0,type="l",lty=2,ylim=c(0,0.2))
#lines(error1,lty=1,col=2) #red yf(x)+random
#lines(error2,lty=1,col=3) #green yf(x)
#if (type=="lovelive"){
#  lines(error3,lty=1,col=4)
#}

