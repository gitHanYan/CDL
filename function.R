###D-trace function
D_trace=function(A,B,lambda,rho){
  #######Input:
  #A,B: covariance matrices of the data matrices
  #lambda: tuning parameter
  #rho: the coefficient of the augmented Lagrangian in D-trace,
  #usually set to 1
  #######Output:
  #D3: result of differential matrix
  maxIter=600
  k=ncol(A)
  D1=solve(diag(diag(A))+diag(k))-solve(diag(diag(B))+diag(k))
  D2=D1
  D3=D1
  
  oldD1=D1+0.01*diag(k); 
  oldD2=D2+0.01*diag(k); 
  oldD3=D3+0.01*diag(k);
  Z=matrix(0,k,k);
  L1=Z;
  L2=Z;
  L3=Z;
  iter=0;
  a1=svd(A)[[1]]#S
  a2=svd(A)[[2]]#U
  b1=svd(B)[[1]]
  b2=svd(B)[[2]]
  temp1=(a1%*%t(b1)+4*rho)
  while (iter<maxIter&&iter%%20!=0|norm((D1-oldD1),"F")/max(1,norm(D1,"F"),norm(oldD1,"F"))>10^-3|
         norm((D2-oldD2),"F")/max(1,norm(D2,"F"),norm(oldD2,"F"))>10^-3|
         norm((D3-oldD3),"F")/max(1,norm(D3,"F"),norm(oldD3,"F"))>10^-3) {
    oldD1=D1
    oldD2=D2
    oldD3=D3
    D1=a2%*%((t(a2)%*%(2*rho*D3+2*rho*D2+A-B+2*L1-2*L3)%*%t(b2))/temp1)%*%b2
    D2=b2%*%((t(b2)%*%(2*rho*D3+2*rho*D1+A-B+2*L3-2*L2)%*%t(a2))/temp1)%*%a2
    D3=sign((rho*D1+rho*D2-L1+L2)/(2*rho))*pmax(abs((rho*D1+rho*D2-L1+L2)/(2*rho))-lambda/(2*rho),0)
    L1=L1+rho*(D3-D1)
    L2=L2+rho*(D2-D3)
    L3=L3+rho*(D1-D2)
    iter=iter+1
  }
  return(D3)
}
CDL=function(A,B,K,lambda,rho){
  #######Input:
  #A, B: the data matrices
  #K: the number of clusters
  #lambda: tuning parameter
  #rho: a tuning parameter in D-trace, usually set to 1
  #######Output:
  #res: result of differential matrix
  res=matrix(0,ncol(A),ncol(A))
  mean_cov=(abs(cov(A))+abs(cov(B)))/2
  Hc_AB=as.dist(1-mean_cov)
  hc_ab=hclust(Hc_AB, method="average",members=NULL)
  cut_AB=cutree(hc_ab,k=K)
  if(1 %in% table(cut_AB)){
    res=D_trace(cov(A),cov(B),lambda,rho)
    return(res)
  }else{
    loc_list=list()
    A_list=list()
    B_list=list()
    for(l in 1:K){
      loc_list[[l]]=which(cut_AB==l)
    }
    for(l in 1:K){
      A_list[[l]]=cov(A[,which(cut_AB==l)])
      B_list[[l]]=cov(B[,which(cut_AB==l)])
    }
    ######ROC curve
    for (l in 1:K) {
      res[loc_list[[l]],loc_list[[l]]]=D_trace(A_list[[l]],B_list[[l]],lambda,rho)
    }
    return(res)
  }
}
