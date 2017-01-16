rm(list=ls(all=TRUE))
library(e1071)
library(ggplot2)
library(quadprog)
####### cannot run if there is no nonsupport vector in an iteration ########
########################## Settings ###########################
#set.seed(11)
nrow = 20000 #put number of rows here to change, data size is nrow/2
dataset = 1 #data set
times = 10 #Number of experiments
type = "radial"  #kernel type
k = 2  #Number added in each iteration
per = 0.005    #Selecting percentage
iter = 20   #the iteration number of adding number

initial.number = 50   #initial set size
cost = 1  #cost variable

###########################################################
K<-function(x,y,kernel="radial",sigma2=2) {
  if(is.vector(x)) {
    if (is.vector(y)) { if (kernel=="linear") { return ( as.matrix(sum(x*y)) ) } else { return ( as.matrix(exp(-sum((x-y)^2)/sigma2)) ) } } else {
      if (kernel=="linear") return ( t(x) %*% t(y) ) else {
        res<-matrix(0,1,nrow(y))
        for ( j in 1:nrow(y) ) { 
          res[1,j]<-exp(-sum((x-y[j,])^2)/sigma2) 
        }
        return (res)
      }
    }
  } else {
    if (is.vector(y)) {
      if (kernel=="linear") return ( x %*% y ) else {
        res<-matrix(0,nrow(x),1)
        for ( i in 1:nrow(x) ) { 
          res[i,1]<-exp(-sum((x[i,]-y)^2)/sigma2) 
        }
        return (res)
      }
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
}

svm.prog.ker=function(xx,yy,lambda=1,type){
  n=nrow(xx)
  p=ncol(xx)
  
  ker.mat=K(xx,xx,type)
  Dmat= (ker.mat+matrix(1,n,n))*(yy%*%t(yy))+diag(1e-5,n)
  dvec=rep(1,n)
  Amat=cbind(diag(1,n),diag(-1,n))
  bvec=c(rep(0,n),rep(-1/lambda,n))
  sol=solve.QP(Dmat,dvec,Amat,bvec,meq=0)
  alpha=sol$solution
  normsq=sum(diag(alpha)%*%Dmat%*%diag(alpha))
  
  list(alpha=alpha,normsq=normsq)
}

delete.row = function(x,position){
  if (is.null(x)) { return(NULL) } else {
    if( dim(x)[1] == 1) {
      return(NULL)
    } else if (dim(x)[1] == 2) {
      return(t(as.matrix(x[-position,])))
    } else if (dim(x)[2] == 1) {
      return(as.matrix(x[-position]))
    } else {
      return(x[-position,])
    }
  }
}

delete.column = function(x,position){
  if (is.null(x)) { return(NULL) } else {
    if( dim(x)[2] == 1) {
      return(NULL)
    } else if (dim(x)[1] == 1) {
      return(t(as.matrix(x[-position])))
    } else {
      return(as.matrix(x[,-position]))
    }
  }
}

initial.solve = function(x,y,samp.set,type){
  m = svm.prog.ker(x,y,lambda=1,type)
  alpha=m$alpha
  tol = 1e-5
  
  inde = samp.set[alpha>cost-tol]    # indices of error vectors; none initially
  indo = samp.set[alpha<=tol]             # indices of other vectors; all initially
  inds = setdiff(setdiff(samp.set,inde),indo)     # indices of support vectors; none initially
  
  l = length(samp.set)
  e=(1:l)[alpha>cost-tol]; o=(1:l)[alpha<=tol]; s=setdiff(setdiff(1:l,e),o)
  
  alpha = alpha[alpha>tol]
  alpha = as.matrix(alpha[alpha<=cost-tol])
  
  Qss = (K(x[s,],x[s,],kernel=type)+1)*(y[s]%*%t(y[s]))
  Qos = (K(x[o,],x[s,],kernel=type)+1)*(y[o]%*%t(y[s]))
  if (length(inde) > 0) {
    Qoe = (K(x[o,],x[e,],kernel=type)+1)*(y[o]%*%t(y[e]))
    Qes = (K(x[e,],x[s,],kernel=type)+1)*(y[e]%*%t(y[s]))
    Qee = (K(x[e,],x[e,],kernel=type)+1)*(y[e]%*%t(y[e]))
    g.e = Qes%*%alpha+rowSums(Qee)-1
    g.o = Qos%*%alpha+rowSums(Qoe)-1
    g.e = -abs(g.e);
    g.o = abs(g.o);
    Q = Qss+diag(1e-5,length(inds))
    R = solve(Q)
    return(list(R=R,Q=Q,g.o=g.o,g.e=g.e,Qos=Qos,Qoe=Qoe,Qes=Qes,Qee=Qee,inds=inds,inde=inde,indo=indo,alpha=alpha))
  } else {
    g.o = Qos%*%alpha-1
    g.o = abs(g.o);
    Q = Qss+diag(1e-5,length(inds))
    R = solve(Q)
    return(list(R=R,Q=Q,g.o=g.o,Qos=Qos,inds=inds,inde=inde,indo=indo,alpha=alpha))
  }
  
}

data.choose = function(x,y,candQ.set,result,k,method){
  inds=result$inds; inde=result$inde; 
  if (length(inds) > 0){
    alpha=result$alpha 
    if (length(inde) > 0){
      cand.fit = ((K(x[candQ.set,],x[inds,],kernel=type)+1)*(y[candQ.set]%*%t(y[inds])))%*%alpha + as.matrix(rowSums((K(x[candQ.set,],x[inde,],kernel=type)+1)*(y[candQ.set]%*%t(y[inde]))))
    } else {
      cand.fit = ((K(x[candQ.set,],x[inds,],kernel=type)+1)*(y[candQ.set]%*%t(y[inds])))%*%alpha
    }
  } else {
    cand.fit = as.matrix(rowSums((K(x[candQ.set,],x[inde,],kernel=type)+1)*(y[candQ.set]%*%t(y[inde]))))
  }
  
  if (method == 1){
    sel.set=candQ.set[order(cand.fit)[1:k]]
    return(sel.set=sel.set)
  } else if (method == 2) {
    rank.set = c(1:length(candQ.set))[cand.fit<1-1e-5]
    if (length(rank.set) > k){
      sel.set=candQ.set[rank.set[order(1-cand.fit[rank.set])[1:k]]]
    } else {
      sel.set=candQ.set[order(cand.fit)[1:k]]
    }
    return(sel.set=sel.set)
  } else if (method == 3) {
    furthest=candQ.set[which.min(cand.fit)]
    t1 = x[candQ.set,]; t2 = matrix(x[furthest,],nrow=1)
    for (i in 1:length(x[furthest,])){
      t1[,i] = t1[,i]-t2[,i]
    }
    sel.set=candQ.set[order(rowSums(t1^2))[2:(k+1)]]
    return(sel.set=sel.set)
  } else if (method == 4) {
    sel.set=candQ.set[order(abs(cand.fit))[1:k]]
    return(sel.set=sel.set)
  } else if (method == 5){
    l = length(candQ.set)
    tuneP = 2; diagK = matrix(0,l)
    for (i in 1:l) {
      diagK[i,] = (K(x[i,],x[i,],kernel=type)+1)*(y[i]%*%t(y[i]))
    }
    sel.set=candQ.set[order(abs(cand.fit/diagK^tuneP))[1:k]]
    return(sel.set=sel.set)
  } else {
    return(sel.set=NULL)
  }
  
}

data.eval = function(x,y,x.test,y.test,result){
  inds=result$inds; inde=result$inde; 
  if (length(inds) > 0) {
    alpha=result$alpha 
    if (length(inde) > 0){
      cand.fit = (K(x.test,x[inds,],kernel=type)+1)%*%(alpha*y[inds]) + as.matrix(rowSums((K(x.test,x[inde,],kernel=type)+1)%*%matrix(y[inde])))
    } else {
      cand.fit = (K(x.test,x[inds,],kernel=type)+1)%*%(alpha*y[inds])
    }
  } else {
    cand.fit = as.matrix(rowSums((K(x.test,x[inde,],kernel=type)+1)%*%matrix(y[inde])))
  }
  return(cand.fit)
}

online.update = function(sel.set,result){
  inds=result$inds; indo=result$indo; inde=result$inde; g.o=result$g.o;
  ls = length(inds); le = length(inde)
  if (le > 0){
    Qoe=result$Qoe; Qee=result$Qee; Qes=result$Qes; g.e=result$g.e
    Qos=NULL; alpha=NULL; R=NULL; Q=NULL;
  }
  if (ls > 0){
    Qos=result$Qos; alpha=result$alpha; R=result$R; Q=result$Q; 
  }
  add.number = length(sel.set)
  
  for (l in 1:add.number){
    ls = length(inds); le = length(inde)
    indc = sel.set[l]
    if (ls > 0){
      if (le > 0){
        g.c = ((K(x[indc,],x[inds,],kernel=type)+1)*(y[indc]%*%t(y[inds])))%*%alpha + as.matrix(rowSums((K(x[indc,],x[inde,],kernel=type)+1)*(y[indc]%*%t(y[inde]))))-1
      } else {
        g.c = ((K(x[indc,],x[inds,],kernel=type)+1)*(y[indc]%*%t(y[inds])))%*%alpha-1
      }
    } else {
      g.c = as.matrix(rowSums((K(x[indc,],x[inde,],kernel=type)+1)*(y[indc]%*%t(y[inde]))))-1
    }
    
    if (g.c > -1e-5){
      if (ls > 0){
        Qsc = (K(x[inds,],x[indc,],kernel=type)+1)*(y[inds]%*%t(y[indc]))
        Qos = rbind(Qos,t(Qsc))
        if (le > 0){
          Qec = (K(x[inde,],x[indc,],kernel=type)+1)*(y[inde]%*%t(y[indc]))
          Qoe = rbind(Qoe,t(Qec))
        }
        g.o = rbind(g.o,abs(g.c))
        indo = c(indo,indc) 
      } else {
        Qsc = NULL; Qos = NULL; 
        Qec = (K(x[inde,],x[indc,],kernel=type)+1)*(y[inde]%*%t(y[indc]))
        Qoe = rbind(Qoe,t(Qec))
        g.o = rbind(g.o,abs(g.c))
        indo = c(indo,indc) 
      }
    } else {
      alphac = 0
      Qcc = (K(x[indc,],x[indc,],kernel=type)+1)*(y[indc]%*%t(y[indc]))
      Qoc = (K(x[indo,],x[indc,],kernel=type)+1)*(y[indo]%*%t(y[indc]))
      if (ls > 0) {
        Qsc = (K(x[inds,],x[indc,],kernel=type)+1)*(y[inds]%*%t(y[indc]))
        b = -R%*%Qsc
        if (le > 0) {
          Qec = (K(x[inde,],x[indc,],kernel=type)+1)*(y[inde]%*%t(y[indc]))
          gamma.e = Qec + Qes%*%b
          gamma.o = Qoc + Qos%*%b
          gamma.c = Qcc + t(Qsc)%*%b
        } else {
          Qec = NULL; gamma.e = NULL
          gamma.o = Qoc + Qos%*%b
          gamma.c = Qcc + t(Qsc)%*%b
        }
      } else {
        Qsc = NULL; b = NULL; 
        Qec = (K(x[inde,],x[indc,],kernel=type)+1)*(y[inde]%*%t(y[indc]))
        gamma.e = Qec
        gamma.o = Qoc
        gamma.c = Qcc
      }
      
      count = 0
      ######################## add #############################
      while (length(indc) > 0) {
        ls = length(inds); le = length(inde);
        alpha1.all = -g.o/gamma.o
        alpha1.p = alpha1.all[gamma.o<0]
        alpha1 = min(c(alpha1.p,cost+100))
        index1 = indo[alpha1.all>=0][which.min(alpha1.p)]
        
        if (le == 0) {
          alpha2 = cost+100
          index2 = 0
        } else {
          alpha2.all = -g.e/gamma.e
          alpha2.p = alpha2.all[gamma.e>0]
          alpha2 = min(c(alpha2.p,cost+100))
          index2 = inde[alpha2.all>=0][which.min(alpha2.p)]
        }
        
        if (ls > 0) {
          alpha3.all = -alpha/b
          alpha3.p = alpha3.all[b<0]
          alpha3 = min(c(alpha3.p,cost+100))
          index3 = inds[alpha3.all>=0][which.min(alpha3.p)]
          
          alpha4.all = (cost-alpha)/b
          alpha4.p = alpha4.all[b>0]
          alpha4 = min(c(alpha4.p,cost+100))
          index4 = inds[alpha4.all>=0][which.min(alpha4.p)]
        } else {
          alpha3 = cost + 100
          alpha4 = cost + 100
        }
        
        alpha5 = -g.c/gamma.c
        alpha6 = cost-alphac
        
        alpha.change = c(alpha1,alpha2,alpha3,alpha4,alpha5,alpha6)
        alpha.p = alpha.change[alpha.change>=0]
        s = alpha.change[alpha.change>=0][which.min(alpha.p)]
        
        ###################### Update ######################
        if (s==alpha1){
          temp = (1:length(indo))[indo==index1]
          index = indo[temp]
          
          g.o = g.o + s*gamma.o
          if (le > 0){
            g.e = g.e + s*gamma.e 
          }
          g.c = g.c + s*gamma.c
          
          Qchange = (K(x[index,],x[index,],kernel=type)+1)*(y[index]%*%t(y[index]))
          Qochange = (K(x[indo,],x[index,],kernel=type)+1)*(y[indo]%*%t(y[index]))
          Qcchange = Qoc[temp,]
          if (ls > 0){
            alpha = alpha + s*b  
            b.change = -R%*%(Qos[temp,])
            gamma.change = Qchange + Qos[temp,]%*%b.change
            R = cbind(rbind(R,rep(0,ls)),rep(0,ls+1)) + c(b.change,1)%*%t(c(b.change,1))/gamma.change[1]    
            Q = cbind(rbind(Q,Qos[temp,]),c(Qos[temp,],Qchange))
          } else {
            gamma.change = Qchange
            R = as.matrix(1/gamma.change[1])
            Q = as.matrix(Qchange)
          }
          Qsc = rbind(Qsc,Qcchange)
          b = -R%*%Qsc
          
          Qos = delete.row((cbind(Qos,Qochange)),temp)
          Qoc = delete.row(Qoc,temp)
          if (le > 0){
            Qes = cbind(Qes,Qoe[temp,])
            Qoe = delete.row(Qoe,temp)
            gamma.e = Qec + Qes%*%b
          }
          gamma.o = Qoc + Qos%*%b
          gamma.c = Qcc + t(Qsc)%*%b
          g.o = delete.row(g.o,temp)
          
          alpha = rbind(alpha,0)
          indo = setdiff(indo,index1)
          inds = c(inds,index1)
          alphac = alphac + s
        } 
        if (s==alpha2){
          temp = (1:le)[inde==index2]
          index = inde[temp]
          
          g.o = g.o + s*gamma.o
          if (le > 0){
            g.e = g.e + s*gamma.e
          }
          g.c = g.c + s*gamma.c
          
          Qchange = (K(x[index,],x[index,],kernel=type)+1)*(y[index]%*%t(y[index]))
          Qochange = (K(x[indo,],x[index,],kernel=type)+1)*(y[indo]%*%t(y[index]))
          Qcchange = Qec[temp,]
          if (ls > 0){
            alpha = alpha + s*b
            b.change = -R%*%(Qes[temp,])
            gamma.change = Qchange + Qes[temp,]%*%b.change
            R = cbind(rbind(R,rep(0,ls)),rep(0,ls+1)) + c(b.change,1)%*%t(c(b.change,1))/gamma.change[1]
            Q = cbind(rbind(Q,Qes[temp,]),c(Qes[temp,],Qchange))
          } else {
            gamma.change = Qchange
            R = as.matrix(1/gamma.change[1])
            Q = as.matrix(Qchange)
          }
          Qsc = rbind(Qsc,Qcchange)
          b = -R%*%Qsc
          
          Qos = cbind(Qos,Qochange)
          Qec = delete.row(Qec,temp)
          Qes = delete.row(cbind(Qes,Qee[temp,]),temp)
          if ( le > 1) {
            Qee = delete.column(delete.row(Qee,temp),temp) 
          } else {
            Qee = NULL
          }
          Qoe = delete.column(Qoe,temp)
          
          g.e = delete.row(g.e,temp)
          if (le > 1) {
            gamma.e = Qec + Qes%*%b
            g.e = -abs(g.e)
          } 
          gamma.o = Qoc + Qos%*%b
          gamma.c = Qcc + t(Qsc)%*%b
          
          alpha = rbind(alpha,cost)
          inde = setdiff(inde,index2)
          inds = c(inds,index2)  
          alphac = alphac + s
        }
        if (s==alpha3){
          temp = (1:ls)[inds==index3]
          index = inds[temp]
          
          g.o = g.o + s*gamma.o
          if (le > 0){
            g.e = g.e + s*gamma.e
          }
          g.c = g.c + s*gamma.c
          
          alpha = alpha + s*b
          alpha = delete.row(alpha,temp)
          k = temp
          for (i in 1:(dim(R)[1])) {
            if (i != k){
              for (j in 1:(dim(R)[2])){
                if (j != k) {
                  R[i,j] = R[i,j] - R[i,k]*R[k,j]/R[k,k]
                }
              }
            }
          }
          R = delete.column(delete.row(R,temp),temp) 
          
          Qoc = rbind(Qoc,Qsc[temp,])
          Qsc = delete.row(Qsc,temp)
          if (ls > 1){
            b = -R%*%Qsc
          }
          if (le > 0){
            Qoe = rbind(Qoe,Qes[,temp])
            Qes = delete.column(Qes,temp)
          }
          Qos = delete.column(rbind(Qos,Q[temp,]),temp)
          Q = delete.column(delete.row(Q,temp),temp) 
          
          if (ls > 1){
            if (length(inde) > 0){
              gamma.e = Qec + Qes%*%b   
            }
            gamma.o = Qoc + Qos%*%b
            gamma.c = Qcc + t(Qsc)%*%b
          } else {
            if (le > 0){
              gamma.e = Qec
            }
            gamma.o = Qoc 
            gamma.c = Qcc 
          }
          
          g.o = rbind(g.o,0)
          inds = setdiff(inds,index)
          indo = c(indo,index)
          alphac = alphac + s
        }
        if (s==alpha4){
          temp = (1:ls)[inds==index4]
          index = inds[temp]
          
          g.o = g.o + s*gamma.o
          if (length(inde) > 0){
            g.e = g.e + s*gamma.e
          }
          g.c = g.c + s*gamma.c
          
          alpha = alpha + s*b
          alpha = delete.row(alpha,temp)
          k = temp
          for (i in 1:(dim(R)[1])) {
            if (i != k){
              for (j in 1:(dim(R)[2])){
                if (j != k) {
                  R[i,j] = R[i,j] - R[i,k]*R[k,j]/R[k,k]
                }
              }
            }
          }
          R = delete.column(delete.row(R,temp),temp) 
          
          Qec = rbind(Qec,Qsc[temp,])
          Qsc = delete.row(Qsc,temp)
          if (ls > 1){
            b = -R%*%Qsc
          }
          if (le > 0 ){
            Qee = cbind(rbind(Qee,Qes[,temp]),c(Qes[,temp],Q[temp,temp]))
            Qoe = cbind(Qoe,Qos[,temp])
            Qes = delete.column(rbind(Qes,Q[temp,]),temp)
          } else {
            Qee = as.matrix(Q[temp,temp])
            Qoe = matrix(Qos[,temp],ncol=1)
            Qes = delete.column(matrix(Q[temp,],nrow=1),temp)
          }
          Qos = delete.column(Qos,temp)
          Q = delete.column(delete.row(Q,temp),temp) 
          
          if (ls > 1){
            gamma.e = Qec + Qes%*%b  
            gamma.o = Qoc + Qos%*%b
            gamma.c = Qcc + t(Qsc)%*%b
          } else {
            gamma.e = Qec
            gamma.o = Qoc
            gamma.c = Qcc
          }
          
          if (le > 0){
            g.e = rbind(g.e,0)
          } else {
            g.e = as.matrix(0)
          }
          
          inds = setdiff(inds,index)
          inde = c(inde,index)
          alphac = alphac + s
        }
        if (s==alpha5){
          Q = rbind(cbind(Q,Qsc),c(Qsc,Qcc))
          Qos = cbind(Qos,Qoc)
          if (le > 0){
            Qes = cbind(Qes,Qec)
            g.e = g.e + s*gamma.e
          }
          g.o = g.o + s*gamma.o
          g.c = g.c + s*gamma.c
          
          alphac = alphac + s
          if (is.null(R)){
            R = as.matrix(1/gamma.c[1])
            alpha = as.matrix(alphac)
          } else {
            R = cbind(rbind(R,rep(0,ls)),rep(0,ls+1)) + c(b,1)%*%t(c(b,1))/gamma.c[1]
            alpha = alpha + s*b
            alpha = rbind(alpha,alphac)
          }
          
          inds = c(inds,indc)
          indc = c()
        }
        if (s==alpha6){
          if (ls > 0){    
            if (le > 0){
              Qoe = cbind(Qoe,as.matrix(Qoc)) 
              Qes = rbind(Qes,t(Qsc))
              Qee = rbind(cbind(Qee,Qec),c(Qec,Qcc))
            } else {
              Qoe = Qoc 
              Qes = t(Qsc)
              Qee = Qcc
            }
          } else {
            Qoe = cbind(Qoe,as.matrix(Qoc)) 
            Qee = rbind(cbind(Qee,Qec),c(Qec,Qcc))
          }
          
          g.o = g.o + s*gamma.o
          g.c = g.c + s*gamma.c
          if (le > 0){
            g.e = g.e + s*gamma.e
            g.e = rbind(g.e,g.c)
          } else {
            g.e = as.matrix(g.c)
          }
          
          alpha = alpha + s*b
          inde = c(inde,indc)
          alphac = alphac + s
          indc = c()
        }
        count = count + 1
      }
    }
  }
  if (length(inde) > 0){
    return(list(R=R,Q=Q,g.o=g.o,g.e=g.e,Qos=Qos,Qoe=Qoe,Qes=Qes,Qee=Qee,inds=inds,inde=inde,indo=indo,alpha=alpha))
  } else{
    return(list(R=R,Q=Q,g.o=g.o,Qos=Qos,inds=inds,inde=inde,indo=indo,alpha=alpha))
  }
}

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
  ptm = proc.time()[1]
}

########################The simulation#####################
seed_num = sample(1:100,times)
error0=0; error1=0; error2=0; error3=0; error4=0; error5=0; error6=0
e0=c(); e1=c(); e2=c(); e3=c(); e4=c(); e5=c(); e6=c()
t = c(0,0,0,0,0,0)
trainsize = nrow(x);

for (tim in 1:times) {
  ############# furthest min yf(x) ############
  set.seed(seed_num[tim])
  samp.set=sample(1:trainsize,initial.number,F)
  cand.set=setdiff(1:trainsize,samp.set)
  test.err1=NULL
  ptm11=0; ptm12=0;
  xx=x[samp.set,]
  yy=y[samp.set]
  result = initial.solve(xx,yy,samp.set,type)

  for(i in 1:iter) {
    ptm_temp <- proc.time()[1]
    result.pred=sign(data.eval(x,y,x.test,y.test,result))
    test.err1=c(test.err1,sum(y.test!=result.pred)/length(y.test))
    ptm11 <- ptm11 + proc.time()[1]-ptm_temp
    
    ptm_temp1 <- proc.time()[1]
    cd = length(cand.set)
    candQ.set=sample(cand.set,cd*per,F)
    sel.set = data.choose(x,y,candQ.set,result,k,1)
    result = online.update(sel.set,result)
    cand.set = setdiff(cand.set,sel.set)
    samp.set = c(samp.set,sel.set)
    ptm12 <- ptm12 + proc.time()[1]-ptm_temp1
  }
  ############# nearest max yf(x)|yf(x)<1 ############
  set.seed(seed_num[tim])
  samp.set=sample(1:trainsize,initial.number,F)
  cand.set=setdiff(1:trainsize,samp.set)
  test.err2=NULL
  ptm21=0; ptm22=0;
  xx=x[samp.set,]
  yy=y[samp.set]
  result = initial.solve(xx,yy,samp.set,type)
  
  for(i in 1:iter) {
    ptm_temp <- proc.time()[1]
    result.pred=sign(data.eval(x,y,x.test,y.test,result))
    test.err2=c(test.err2,sum(y.test!=result.pred)/length(y.test))
    ptm21 <- ptm21 + proc.time()[1]-ptm_temp
    
    ptm_temp1 <- proc.time()[1]
    cd = length(cand.set)
    candQ.set=sample(cand.set,cd*per,F)
    sel.set = data.choose(x,y,candQ.set,result,k,2)
    result = online.update(sel.set,result)
    cand.set = setdiff(cand.set,sel.set)
    samp.set = c(samp.set,sel.set)
    ptm22 <- ptm22 + proc.time()[1]-ptm_temp1
  }
  ############# KNN furthest min yf(x)|near to min yf(x) ############
  set.seed(seed_num[tim])
  samp.set=sample(1:trainsize,initial.number,F)
  cand.set=setdiff(1:trainsize,samp.set)
  test.err3=NULL
  ptm31=0; ptm32=0;
  xx=x[samp.set,]
  yy=y[samp.set]
  result = initial.solve(xx,yy,samp.set,type)
  
  for(i in 1:iter) {
    ptm_temp <- proc.time()[1]
    result.pred=sign(data.eval(x,y,x.test,y.test,result))
    test.err3=c(test.err3,sum(y.test!=result.pred)/length(y.test))
    ptm31 <- ptm31 + proc.time()[1]-ptm_temp
    
    ptm_temp1 <- proc.time()[1]
    cd = length(cand.set)
    candQ.set=sample(cand.set,cd*per,F)
    sel.set = data.choose(x,y,candQ.set,result,k,3)
    result = online.update(sel.set,result)
    cand.set = setdiff(cand.set,sel.set)
    samp.set = c(samp.set,sel.set)
    ptm32 <- ptm32 + proc.time()[1]-ptm_temp1
  }
  ############# |yf(x)| ###################
  set.seed(seed_num[tim])
  samp.set=sample(1:trainsize,initial.number,F)
  cand.set=setdiff(1:trainsize,samp.set)
  test.err4=NULL
  ptm41=0; ptm42=0;
  xx=x[samp.set,]
  yy=y[samp.set]
  result = initial.solve(xx,yy,samp.set,type)
  
  for(i in 1:iter) {
    ptm_temp <- proc.time()[1]
    result.pred=sign(data.eval(x,y,x.test,y.test,result))
    test.err4=c(test.err4,sum(y.test!=result.pred)/length(y.test))
    ptm41 <- ptm41 + proc.time()[1]-ptm_temp
    
    ptm_temp1 <- proc.time()[1]
    cd = length(cand.set)
    candQ.set=sample(cand.set,cd*per,F)
    sel.set = data.choose(x,y,candQ.set,result,k,4)
    result = online.update(sel.set,result)
    cand.set = setdiff(cand.set,sel.set)
    samp.set = c(samp.set,sel.set)
    ptm42 <- ptm42 + proc.time()[1]-ptm_temp1
  }
  ############# |yf(x)|/K ################
  set.seed(seed_num[tim])
  samp.set=sample(1:trainsize,initial.number,F)
  cand.set=setdiff(1:trainsize,samp.set)
  test.err5=NULL
  ptm51=0; ptm52=0;
  xx=x[samp.set,]
  yy=y[samp.set]
  result = initial.solve(xx,yy,samp.set,type)

  for(i in 1:iter) {
    ptm_temp <- proc.time()[1]
    result.pred=sign(data.eval(x,y,x.test,y.test,result))
    test.err5=c(test.err5,sum(y.test!=result.pred)/length(y.test))
    ptm51 <- ptm51 + proc.time()[1]-ptm_temp
    
    ptm_temp1 <- proc.time()[1]
    cd = length(cand.set)
    candQ.set=sample(cand.set,cd*per,F)
    sel.set = data.choose(x,y,candQ.set,result,k,5)
    result = online.update(sel.set,result)
    cand.set = setdiff(cand.set,sel.set)
    samp.set = c(samp.set,sel.set)
    ptm52 <- ptm52 + proc.time()[1]-ptm_temp1
  } 
  ############# random ###################
  set.seed(seed_num[tim])
  samp.set=sample(1:trainsize,initial.number,F)
  cand.set=setdiff(1:trainsize,samp.set)
  ptm01=0; ptm02=0;
  test.err0=NULL
  xx=x[samp.set,]
  yy=y[samp.set]
  result = initial.solve(xx,yy,samp.set,type)
  
  for(i in 1:iter) {
    ptm_temp <- proc.time()[1]
    result.pred=sign(data.eval(x,y,x.test,y.test,result))
    test.err0=c(test.err0,sum(y.test!=result.pred)/length(y.test))
    ptm01 <- ptm01 + proc.time()[1]-ptm_temp
    
    ptm_temp1 <- proc.time()[1]
    sel.set=sample(cand.set,k,F)
    result = online.update(sel.set,result)
    cand.set = setdiff(cand.set,sel.set)
    samp.set = c(samp.set,sel.set)
    ptm02 <- ptm02 + proc.time()[1]-ptm_temp1
  }
  
  error0 = error0+test.err0; error1 = error1+test.err1; error2 = error2+test.err2
  error3 = error3+test.err3; error4 = error4+test.err4; error5 = error5+test.err5
  e0 = rbind(e0,test.err0); e1 = rbind(e1,test.err1); e2 = rbind(e2,test.err2)
  e3 = rbind(e3,test.err3); e4 = rbind(e4,test.err4); e5 = rbind(e5,test.err5)
  
  t = t + c(ptm02,ptm12,ptm22,ptm32,ptm42,ptm52)
}

t = t/times
############### direct solve ##################
m=svm(x,as.factor(y),cost=1,kernel=type)
m.pred=predict(m,x.test,decision.values=T)
base = sum(y.test!=m.pred)/length(y.test)
base = rep(base,iter)
############### Analyse error ################
se0=c(); se1=c(); se2=c(); se3=c(); se4=c(); se5=c(); be=rep(0,iter)
for (i in 1:iter){
  se0 = c(se0,sd(e0[,i])/sqrt(times))
  se1 = c(se1,sd(e1[,i])/sqrt(times))
  se2 = c(se2,sd(e2[,i])/sqrt(times))
  se3 = c(se3,sd(e3[,i])/sqrt(times))
  se4 = c(se4,sd(e4[,i])/sqrt(times))
  se5 = c(se5,sd(e5[,i])/sqrt(times))
}

error0 = error0/times; error1 = error1/times; error2 = error2/times
error3 = error3/times; error4 = error4/times; error5 = error5/times

################## Plot #######################
l = iter
iteration = c(rep(1:l),rep(1:l),rep(1:l),rep(1:l),rep(1:l),rep(1:l),rep(1:l))
error = c(error0,error1,error2,error3,error4,error5,base)
color_label = c(rep("random",l),rep("furthest",l),rep("nesrest",l),rep("furthest KNN",l),rep("|yf(x)|",l),rep("yf(x)/K",l),rep("base",l))
se = c(se0,se1,se2,se3,se4,se5,be)
e = data.frame(iteration,error,color_label,se)
pd <- position_dodge(0.2)

#ggplot(e, aes(x=iteration, y=error, colour=color_label)) + 
#  geom_errorbar(aes(ymin=error-se, ymax=error+se), group = color_label, width=.2, position=pd) +
#  geom_line(position=pd) +
#  geom_point(position=pd)

ggplot(e, aes(x=iteration, y=error, colour=color_label)) + 
  geom_errorbar(aes(ymin=error-se, ymax=error+se), width=.2, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd)
