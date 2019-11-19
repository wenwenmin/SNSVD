source('SNSVD_L0SVD_ANSVD_SyntheticTest.R')
# -------------------------------------------
p = 200
q = 100
t = 50 
lambda = seq(0.02,0.06,by = 0.005) 
p11 = 0.3; p12 = 0.1
A1 = get_sim_prior_Net(p,t,p11,p12) # Generate simple network for row variables.
A2 = get_sim_prior_Net(q,t,p11,p12) # Generate simple network for col variables.

SNR = rep(0,length(lambda))
# ----------------------------------------------
# ----------------------------------------------
Res50 = list()
for(kk in 1:50){
  print(kk)
  
  Res = list()
  specificity.u = Sensitivity.u = accuracy.u = matrix(0,nrow = 4, ncol = length(lambda))
  specificity.v = Sensitivity.v = accuracy.v = matrix(0,nrow = 4, ncol = length(lambda))
  singular.value = matrix(0,nrow = 4, ncol = length(lambda))
  
  for(j in 1:length(lambda)){
    snr.lambda = lambda[j]
    u <- matrix(c(rep(1,50), rep(0,p-t)), ncol=1);u <- u/norm(u,'E')
    v <- matrix(c(rep(1,50), rep(0,q-t)), ncol=1);v <- v/norm(v,'E')
    
    set.seed(100*kk)
    X <- u%*%t(v) + snr.lambda*matrix(rnorm(p*q),ncol=q)
    u.true = u; v.true = v; 
    
    SNR[j] = norm(u%*%t(v),"F")^2/norm(snr.lambda*matrix(rnorm(p*q),ncol=q),"F")^2
    # ----------------------------------------------
    # ----------------------------------------------
    ku = t; kv = t
    delta1 = delta2 = 0.5
    
    out1 = SNSVD(X, A1, delta1, ku, kv, niter=200)
    out2 = L0_svd(X, ku, kv, niter=200)
    out3 = ALasso_svd(X, ku, kv, niter=200)
    out4 = SCAD_svd(X, ku, kv, niter=200)
    # ----------------------------------------------
    # ----------------------------------------------
    # u
    Sensitivity.u[1,j] = length(which(out1$u[1:t,1]!=0))/t
    specificity.u[1,j] = length(which(out1$u[(t+1):p,1]==0))/(p-t)
    accuracy.u[1,j] = (length(which(out1$u[1:t,1]!=0)) + length(which(out1$u[(t+1):p,1]==0)))/p
    
    Sensitivity.u[2,j] = length(which(out2$u[1:t,1]!=0))/t
    specificity.u[2,j] = length(which(out2$u[(t+1):p,1]==0))/(p-t)
    accuracy.u[2,j] = (length(which(out2$u[1:t,1]!=0)) + length(which(out2$u[(t+1):p,1]==0)))/p
    
    Sensitivity.u[3,j] = length(which(out3$u[1:t,1]!=0))/t
    specificity.u[3,j] = length(which(out3$u[(t+1):p,1]==0))/(p-t)
    accuracy.u[3,j] = (length(which(out3$u[1:t,1]!=0)) + length(which(out3$u[(t+1):p,1]==0)))/p
    
    Sensitivity.u[4,j] = length(which(out4$u[1:t,1]!=0))/t
    specificity.u[4,j] = length(which(out4$u[(t+1):p,1]==0))/(p-t)
    accuracy.u[4,j] = (length(which(out4$u[1:t,1]!=0)) + length(which(out4$u[(t+1):p,1]==0)))/p
    
    # v
    Sensitivity.v[1,j] = length(which(out1$v[1:t,1]!=0))/t
    specificity.v[1,j] = length(which(out1$v[(t+1):q,1]==0))/(q-t)
    accuracy.v[1,j] = (length(which(out1$v[1:t,1]!=0)) + length(which(out1$v[(t+1):q,1]==0)))/q
    
    Sensitivity.v[2,j] = length(which(out2$v[1:t,1]!=0))/t
    specificity.v[2,j] = length(which(out2$v[(t+1):q,1]==0))/(q-t)
    accuracy.v[2,j] = (length(which(out2$v[1:t,1]!=0)) + length(which(out2$v[(t+1):q,1]==0)))/q
    
    Sensitivity.v[3,j] = length(which(out3$v[1:t,1]!=0))/t
    specificity.v[3,j] = length(which(out3$v[(t+1):q,1]==0))/(q-t)
    accuracy.v[3,j] = (length(which(out3$v[1:t,1]!=0)) + length(which(out3$v[(t+1):q,1]==0)))/q
    
    Sensitivity.v[4,j] = length(which(out4$v[1:t,1]!=0))/t
    specificity.v[4,j] = length(which(out4$v[(t+1):q,1]==0))/(q-t)
    accuracy.v[4,j] = (length(which(out4$v[1:t,1]!=0)) + length(which(out4$v[(t+1):q,1]==0)))/q
    
    singular.value[1,j] = out1$d
    singular.value[2,j] = out2$d
    singular.value[3,j] = out3$d
    singular.value[4,j] = out4$d
  }
  
  Res$Sensitivity.u = Sensitivity.u
  Res$specificity.u = specificity.u
  Res$accuracy.u = accuracy.u
  
  Res$Sensitivity.v = Sensitivity.v
  Res$specificity.v = specificity.v
  Res$accuracy.v = accuracy.v
  
  Res50[[kk]] = Res
}