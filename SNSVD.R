# Sparse Network Constrained SVD (SNSVD)
SNSVD = function(X, A1, delta1, ku, kv, niter=1000, err=0.0001, Num.init = 5){
  # X is a matrix with p x n. There are n samples and p features (genes)
  # A1: adjacency matrix for row variables (genes)
  
  # Compute normalized Laplacian Matrix 
  W1 = laplacian.matrix(A1)
  d.opt = -1000
  for(ii in 1:Num.init){
    set.seed(ii*1000)
    v0 = matrix(rnorm(ncol(X),0,1),ncol=1);v0 = v0/norm(v0,'E')
    # set.seed(200)
    u0 = matrix(rnorm(nrow(X),0,1),ncol=1);u0 = u0/norm(u0,'E')
    
    for (i in 1:niter){
 
      # update u vector
      # u0 = abs(u0)
      u = Nproject(X%*%v0,   u0, ku, W1, delta1) 
      
      # update v vector 
      v = l0project(t(X)%*%u, kv) 
      
      # Algorithm termination condition norm(matrix(v),'E')
      if ((norm(u - u0,'E')<= err)&(norm(v - v0,'E')<= err)){break}
      else {u0 = u; v0 = v}
    }
    d = t(u)%*%X%*%v
    if(d>d.opt){
      d.opt = d
      u.opt = u
      v.opt = v
    }
  }
  return (list(u=u.opt, v=v.opt, d=d.opt))
}

# L0 constrained SVD (L0_svd)
L0_svd = function(X, ku, kv, niter=1000, err=0.0001, Num.init = 5){
  # X: input matrix
  d.opt = 0
  for(ii in 1:Num.init){
    
    set.seed(ii*100)
    v0 = matrix(rnorm(ncol(X),0,1),ncol=1);v0 = v0/norm(v0,'E')
    
    set.seed(ii*200)
    u0 = matrix(rnorm(nrow(X),0,1),ncol=1);u0 = u0/norm(u0,'E')
    
    # Iterative algorithm to solve u and v values
    for (i in 1:niter){
      
      u = l0project(X%*%v0, ku) 
      v = l0project(t(X)%*%u, kv) 
      
      # Algorithm termination condition norm(matrix(v),'E')
      if ((norm(u - u0,'E')<= err)&(norm(v - v0,'E')<= err)){break}
      else {
        u0 = u;v0 = v}
    }
    d = t(u)%*%X%*%v
    if(d>d.opt){
      d.opt = d
      u.opt = u
      v.opt = v
    }
  }
  return (list(u=u.opt, v=v.opt, d=d.opt))
}

ALasso_svd = function(X, ku, kv, niter=1000, err=0.0001, Num.init = 5){
  # X: input matrix
  d.opt = -1000
  for(ii in 1:Num.init){
    
    set.seed(ii*100)
    v0 = matrix(rnorm(ncol(X),0,1),ncol=1);v0 = v0/norm(v0,'E')
    # set.seed(ii*200)
    u0 = matrix(rnorm(nrow(X),0,1),ncol=1);u0 = u0/norm(u0,'E')
    
    # Iterative algorithm to solve u and v values
    for (i in 1:niter){
      
      u = ALasso.project(X%*%v0, ku) 
      v = ALasso.project(t(X)%*%u, kv) 
      
      # Algorithm termination condition norm(matrix(v),'E')
      if ((norm(u - u0,'E')<= err)&(norm(v - v0,'E')<= err)){break}
      else {
        u0 = u;v0 = v}
    }
    d = t(u)%*%X%*%v
    if(d>d.opt){
      d.opt = d
      u.opt = u
      v.opt = v
    }
  }
  return (list(u=u.opt, v=v.opt, d=d.opt))
}

SCAD_svd = function(X, ku, kv, niter=1000, err=0.0001, Num.init = 5){
  # X: input matrix
  d.opt = -1000
  for(ii in 1:Num.init){
    
    set.seed(ii*100)
    v0 = matrix(rnorm(ncol(X),0,1),ncol=1);v0 = v0/norm(v0,'E')
    
    set.seed(ii*200)
    u0 = matrix(rnorm(nrow(X),0,1),ncol=1);u0 = u0/norm(u0,'E')
    
    # Iterative algorithm to solve u and v values
    for (i in 1:niter){
      
      u = SCAD.project(X%*%v0, ku) 
      v = SCAD.project(t(X)%*%u, kv) 
      
      # Algorithm termination condition norm(matrix(v),'E')
      if ((norm(u - u0,'E')<= err)&(norm(v - v0,'E')<= err)){break}
      else {
        u0 = u;v0 = v}
    }
    d = t(u)%*%X%*%v
    if(d>d.opt){
      d.opt = d
      u.opt = u
      v.opt = v
    }
  }
  return (list(u=u.opt, v=v.opt, d=d.opt))
}

# -----------------------------------------------
# network constrained SVD project function
Nproject = function(z, u0, k, W, delta){
  # 2016-5-6  
  z = z + delta*W%*%u0
  u = select2(abs(z),k)
  
  if(sum(u^2)==0){return(rep(0,length(u)))}
  else{
    u = u/sqrt(sum(u^2))
    u = sign(z)*u
    return(u)} 
}

L1Nproject = function(z, u0, k, W, delta){ 
  # 2016-5-6
  z = z + delta*W%*%u0
  lambda = sort(abs(z),decreasing=T)[k+1]
  u = abs(z)-lambda
  u[u<0] <- 0
  
  if(sum(u^2)==0){return(rep(0,length(u)))}
  else{
    u = u/sqrt(sum(u^2))
    u = sign(z)*u
    return(u)} 
}
# L0 constrained SVD project function
l0project = function(z, k){  
  absz = abs(z);
  u = select2(absz,k)
  
  if(sum(u^2)==0){return(rep(0,length(u)))}
  else{
    u = u/sqrt(sum(u^2))
    u = sign(z)*u
    return(u)} 
}
# An auxiliary function
select2 = function(x, k){
  if(k>=length(x)) return(x)
  x[-order(x,decreasing=T)[1:k]] = 0
  return(x)
}

ALasso.project = function(z, k){  
  absz = abs(z); 
  sorz = matrix(sort(absz,decreasing = TRUE))
  delta = sorz[k+1]
  u = sign(z)*(abs(z)>=delta)*(abs(z)-delta)
  
  if(sum(u^2)==0){return(rep(0,length(u)))}
  else{
    u = u/sqrt(sum(u^2))
    return(u)} 
}

SCAD.project = function(z, k){ 
  absz = abs(z); 
  sorz = matrix(sort(absz,decreasing = TRUE))
  delta = sorz[k+1]  
  a = 3.7 # a: default choice for SCAD penalty
  
  u = sign(z)*(abs(z)>=delta)*(abs(z)-delta)*(abs(z)<=2*delta)+
    ((a-1)*z-sign(z)*a*delta)/(a-2)*(2*delta<abs(z))*(abs(z)<=a*delta)+z*(abs(z)>a*delta)
  
  if(sum(u^2)==0){return(rep(0,length(u)))}
  else{
    u = u/sqrt(sum(u^2))
    return(u)} 
}
# ----------------------------------------------
# ----------------------------------------------  
# get_sim_prior_Net(10,5,0.8,0.2)
get_sim_prior_Net = function(n,t,p11,p12){
  A = matrix(0,nrow=n,ncol=n) 
  for(i in 1:n){
    for(j in 1:n){
      if(i>j){
        set.seed(10*i+8*j)
        if(i<t&j<t){
          if(runif(1)<p11) A[i,j] = 1}
        else{
          if(runif(1)<p12) A[i,j] = 1}  
      }
    }
  }
  A = A +t(A)
  diag(A)=0
  return(A)
}
#---------------------------------------------------
#---------------------------------------------------
replace.NAs = function(X){
  Xnew = X
  na.id = is.na(X)
  if(sum(na.id)>0) Xnew[na.id] = mean(X[!na.id])
  return(Xnew)}
#---------------------------------------------------
#---------------------------------------------------
laplacian.matrix = function(A){
  diag(A) <- 0
  A_deg <- apply(abs(A),1,sum)
  p <- ncol(A)
  L <- matrix(0,p,p)
  A_nonzero <- which(A_deg!=0)
  for (I in A_nonzero){
    for (J in A_nonzero){
      L[I,J] <- A[I,J]/sqrt(A_deg[I]*A_deg[J])
    }
  }
  return(L)
}
#---------------------------------------------------
#---------------------------------------------------