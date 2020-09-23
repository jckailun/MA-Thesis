##### Newton's method for cohorts #####
theta.1 = rep(0, 48+2)
theta.2 = rep(0, 48+7)
### FYR, we always resort to these outcome dummies 
D11 = ifelse(d1$acgrd_9 <= 11, 1, 0) # dummy variable indicating that they did not transition at all from 1
D21 = ifelse(d1$acgrd_9 == 12, 1, 0) # dummy variable indicating that they transitioned from 1 - 2
D31 = ifelse((d1$acgrd_9 >= 13 & d1$acgrd_9 <= 15), 1, 0) # dummy variable indicating that they transitioned from 2 - 3
D41 = 1 - D11 - D21 - D31 # dummy variable indicating that they transitioned from 3 - 4

D12 = ifelse(d2$acgrd_10 <= 11, 1, 0) # dummy variable indicating that they did not transition at all from 1
D22 = ifelse(d2$acgrd_10 == 12, 1, 0) # dummy variable indicating that they transitioned from 1 - 2
D32 = ifelse((d2$acgrd_10 >= 13 & d2$acgrd_10 <= 15), 1, 0) # dummy variable indicating that they transitioned from 2 - 3
D42 = 1 - D12 - D22 - D32 # dummy variable indicating that they transitioned from 3 - 4

dummy.1 = cbind(D11, D21, D31, D41)
dummy.2 = cbind(D12, D22, D32, D42)

dummy = list(dummy.1, dummy.2)
###

ll = function(theta, df, k = 16, cohort, dummy){ # validated
  
  G = length(df)
  b1 = theta[1:k]
  b2 = theta[(k+1):(2*k)]
  b3 = theta[(2*k+1):(3*k)]
  a = theta[-(1:(3*k))]
  l = 0
  
  D = dummy[[cohort]]
  #else {
  #  if (G == 2) { D = dummy[[1]] } else {D = dummy[[2]]} }
  
  D1 = D[,1]; D2 = D[,2]; D3 = D[,3]; D4 = D[,4]
  
  for (j in 1:G){
    
    d = as.matrix(df[[j]])
    aj = a[j]
    n = nrow(d)
    
    if (G == 1){ range = 1:n } else {
      if (cohort == 1){ range = idx.d1[[j]] }
      if (cohort == 2){ range = idx.d2[[j]] }
    }
    D11 = D1[range]; D21 = D2[range]; D31 = D3[range]; D41 = D4[range]
    
    for (i in 1:n){
      
      x_i = as.vector(d[i,])
      eta1 = as.numeric(x_i %*% b1 + aj)
      eta2 = as.numeric(x_i %*% b2 + aj)
      eta3 = as.numeric(x_i %*% b3 + aj)
      
      t1 = -log(1 + exp(eta1))*D11[i]
      t2 = (eta1 - log(1 + exp(eta1)) - log(1 + exp(eta2))) * D21[i]
      t3 = (eta1 + eta2 - log(1 + exp(eta1)) - log(1 + exp(eta2)) - log(1 + exp(eta3))) * D31[i]
      t4 = (eta1 + eta2 + eta3 - log(1 + exp(eta1)) - log(1 + exp(eta2)) - log(1 + exp(eta3))) * D41[i]
      
      l = l + (t1 + t2 + t3 + t4)

    }
  }
  l
}
ll(rep(0,49), k=16, list(X1.c), cohort=1, dummy)
### Score Function ###
dldb1 = function(theta, df, k = 16, cohort, dummy){ # validated
  
  G = length(df)
  b1 = theta[1:k]
  b2 = theta[(k+1):(2*k)]
  b3 = theta[(2*k+1):(3*k)]
  a = theta[-(1:(3*k))]
  B1 = 0
  
  D = dummy[[cohort]]
  
  D1 = D[,1]; D2 = D[,2]; D3 = D[,3]; D4 = D[,4]
  
  for (j in 1:G){
    
    d = as.matrix(df[[j]])
    aj = a[j]
    n = nrow(d)
    if (G == 1){ range = 1:n } else {
      if (cohort == 1){ range = idx.d1[[j]] }
      if (cohort == 2){ range = idx.d2[[j]] }
    }
    D11 = D1[range]; D21 = D2[range]; D31 = D3[range]; D41 = D4[range]

    for (i in 1:n){
      x_i = as.vector(d[i,])
      eta1 = as.numeric(x_i %*% b1 + aj)
      #eta2 = as.numeric(x_i %*% b2 + aj)
      #eta3 = as.numeric(x_i %*% b3 + aj)
      
      t1 = -plogis(eta1) * x_i * D11[i]
      t2 = (1-plogis(eta1)) * x_i * D21[i]
      t3 = (1-plogis(eta1)) * x_i * D31[i]
      t4 = (1-plogis(eta1)) * x_i * D41[i]
      
      B1 = B1 + (t1 + t2 + t3 + t4)
    }
  }
  
  as.vector(B1)
  
}
dldb1(rep(0,49), df=list(X1.c), k=16, cohort=1, dummy)

dldb2 = function(theta, df, k=16, cohort, dummy){ 
  
  G = length(df)
  b1 = theta[1:k]
  b2 = theta[(k+1):(2*k)]
  b3 = theta[(2*k+1):(3*k)]
  a = theta[-(1:(3*k))]
  B2 = 0
  
  D = dummy[[cohort]]
  
  D1 = D[,1]; D2 = D[,2]; D3 = D[,3]; D4 = D[,4]
  
  for (j in 1:G){
    
    d = as.matrix(df[[j]])
    aj = a[j]
    n = nrow(d)
    
    if (G == 1){ range = 1:n } else {
      if (cohort == 1){ range = idx.d1[[j]] }
      if (cohort == 2){ range = idx.d2[[j]] }
    }
    
    D11 = D1[range]; D21 = D2[range]; D31 = D3[range]; D41 = D4[range]
      
    for (i in 1:n){
      x_i = as.vector(d[i,])
      #eta1 = as.numeric(x_i %*% b1 + aj)
      eta2 = as.numeric(x_i %*% b2 + aj)
      #eta3 = as.numeric(x_i %*% b3 + aj)
        
      t2 = -plogis(eta2) * x_i * D21[i]
      t3 = (1-plogis(eta2)) * x_i * D31[i]
      t4 = (1-plogis(eta2)) * x_i * D41[i]
        
      B2 = B2 + (t2 + t3 + t4)
    }
  }
    
  as.vector(B2)
    
}
dldb2(rep(0,49), df=list(X1.c), k=16, cohort=1, dummy)

dldb3 = function(theta, df, k=16, cohort, dummy){ # validated
  
  G = length(df)
  b1 = theta[1:k]
  b2 = theta[(k+1):(2*k)]
  b3 = theta[(2*k+1):(3*k)]
  a = theta[-(1:(3*k))]
  B3 = 0
  
  D = dummy[[cohort]]
  
  D1 = D[,1]; D2 = D[,2]; D3 = D[,3]; D4 = D[,4]
  
  for (j in 1:G){
    
    d = as.matrix(df[[j]])
    aj = a[j]
    n = nrow(d)
    
    if (G == 1){ range = 1:n } else {
      if (cohort == 1){ range = idx.d1[[j]] }
      if (cohort == 2){ range = idx.d2[[j]] }
    }
    D11 = D1[range]; D21 = D2[range]; D31 = D3[range]; D41 = D4[range]
    
    for (i in 1:n){
      x_i = as.vector(d[i,])
      #eta1 = as.numeric(x_i %*% b1 + aj)
      #eta2 = as.numeric(x_i %*% b2 + aj)
      eta3 = as.numeric(x_i %*% b3 + aj)
      
      t3 = -plogis(eta3) * x_i * D31[i]
      t4 = (1-plogis(eta3)) * x_i * D41[i]
      
      B3 = B3 + (t3 + t4)
    }
  }
  
  as.vector(B3)
  
}
dldb3(rep(0,49), df=list(X1.c), k=16, cohort=1, dummy)

dldag = function(theta, df, k = 16, g, cohort, dummy){ # validated
  # g controls groups
  G = length(df)
  b1 = theta[1:k]
  b2 = theta[(k+1):(2*k)]
  b3 = theta[(2*k+1):(3*k)]
  aj = theta[3*k+g]
  Ag = 0
  
  D = dummy[[cohort]]
  
  D1 = D[,1]; D2 = D[,2]; D3 = D[,3]; D4 = D[,4]
  
  d = as.matrix(df[[g]])
  n = nrow(d)
  
  if (G == 1){ range = 1:n } else {
    if (cohort == 1){ range = idx.d1[[g]] }
    if (cohort == 2){ range = idx.d2[[g]] }
  }

  D11 = D1[range]; D21 = D2[range]; D31 = D3[range]; D41 = D4[range]

  
  for (i in 1:n){
    x_i = as.vector(d[i,])
    eta1 = as.numeric(x_i %*% b1 + aj)
    eta2 = as.numeric(x_i %*% b2 + aj)
    eta3 = as.numeric(x_i %*% b3 + aj)
    
    t1 = - plogis(eta1) * D11[i]
    t2 = (1 - plogis(eta1) - plogis(eta2))  * D21[i]
    t3 = (2 - plogis(eta1) - plogis(eta2) - plogis(eta3)) * D31[i]
    t4 = (3 - plogis(eta1) - plogis(eta2) - plogis(eta3)) * D41[i]
    
    Ag = Ag + (t1 + t2 + t3 + t4)
  }
  Ag
}
dldag(rep(0,49), df=list(X1.c), k=16, g=1, cohort=1, dummy)

dldag(rep(0,50), df=G1, k=16, g=1, cohort=1, dummy)
dldag(rep(0,51), df=G2, k=16, g=1, cohort=2, dummy)

score = function(theta, df, k=16, cohort, dummy){ # validated
  
  G = length(df)
  sa = NULL
  
  sb1 = dldb1(theta, df, k = 16, cohort, dummy)
  sb2 = dldb2(theta, df, k = 16, cohort, dummy)
  sb3 = dldb3(theta, df, k = 16, cohort, dummy)
  
  for (j in 1:G){
    sa = c(sa, dldag(theta, df, k = 16, g = j, cohort, dummy))
  }
  
  c(sb1, sb2, sb3, sa)
  
}
score(rep(0,50), list(X1.c), k=16, cohort=1, dummy)
### End: Score Function ###

### Hessian Matrix ###
d2ldb12 = function(theta, df, k = 16, dummy){ # validated
  
  G = length(df)
  b1 = theta[1:k]
  b2 = theta[(k+1):(2*k)]
  b3 = theta[(2*k+1):(3*k)]
  a = theta[-(1:(3*k))]
  BB1 = 0
  
  #D = dummy[[cohort]]
  
  #D1 = D[,1]; D2 = D[,2]; D3 = D[,3]; D4 = D[,4]
  
  for (j in 1:G){
    
    d = as.matrix(df[[j]])
    ag = a[j]
    n = nrow(d)
    
#    if (G == 1){ range = 1:n } else {
#      if (cohort == 1){ range = idx.d1[[j]] }
#      if (cohort == 2){ range = idx.d2[[j]] }
#    }
#    D11 = D1[range]; D21 = D2[range]; D31 = D3[range]; D41 = D4[range]
    
    for (i in 1:n){
      
      x_i = as.vector(d[i,])
      eta1 = as.numeric(x_i %*% b1 + ag)
      #eta2 = as.numeric(x_i %*% b2 + aj)
      #eta3 = as.numeric(x_i %*% b3 + aj)
      
      t1 = plogis(eta1) * (1-plogis(eta1)) * (x_i%o%x_i)
      BB1 = BB1 + t1
      
    }
    
  }
  
  -BB1

}
d2ldb12(rep(0,51), k=16, G2, dummy)

d2ldb22 = function(theta, df, k = 16, cohort, dummy){ # validated
  
  G = length(df)
  b1 = theta[1:k]
  b2 = theta[(k+1):(2*k)]
  b3 = theta[(2*k+1):(3*k)]
  a = theta[-(1:(3*k))]
  BB2 = 0
  
  D = dummy[[cohort]]
  
  D1 = D[,1]; D2 = D[,2]; D3 = D[,3]; D4 = D[,4]
  
  for (j in 1:G){
    
    d = as.matrix(df[[j]])
    ag = a[j]
    n = nrow(d)
    
    if (G == 1){ range = 1:n } else {
      if (cohort == 1){ range = idx.d1[[j]] }
      if (cohort == 2){ range = idx.d2[[j]] }
    }
    D11 = D1[range]; D21 = D2[range]; D31 = D3[range]; D41 = D4[range]
    
    for (i in 1:n){
      
      x_i = as.vector(d[i,])
      #eta1 = as.numeric(x_i %*% b1 + aj)
      eta2 = as.numeric(x_i %*% b2 + ag)
      #eta3 = as.numeric(x_i %*% b3 + aj)
      
      t1 = (1 - D11[i]) * plogis(eta2) * (1-plogis(eta2)) * (x_i%o%x_i)
      BB2 = BB2 + t1
      
    }
    
  }
  
  -BB2
  
}
d2ldb22(rep(0,50), k=16, G1, cohort = 1, dummy)

d2ldb32 = function(theta, df, k = 16, cohort, dummy){ # validated
  
  G = length(df)
  b1 = theta[1:k]
  b2 = theta[(k+1):(2*k)]
  b3 = theta[(2*k+1):(3*k)]
  a = theta[-(1:(3*k))]
  BB3 = 0
  
  D = dummy[[cohort]]
  
  D1 = D[,1]; D2 = D[,2]; D3 = D[,3]; D4 = D[,4]
  
  for (j in 1:G){
    
    d = as.matrix(df[[j]])
    ag = a[j]
    n = nrow(d)
    
    if (G == 1){ range = 1:n } else {
      if (cohort == 1){ range = idx.d1[[j]] }
      if (cohort == 2){ range = idx.d2[[j]] }
    }
    D11 = D1[range]; D21 = D2[range]; D31 = D3[range]; D41 = D4[range]
    
    for (i in 1:n){
      
      x_i = as.vector(d[i,])
      #eta1 = as.numeric(x_i %*% b1 + aj)
      #eta2 = as.numeric(x_i %*% b2 + aj)
      eta3 = as.numeric(x_i %*% b3 + ag)
      
      t1 = (D31[i]+D41[i]) * plogis(eta3) * (1-plogis(eta3)) * (x_i%o%x_i)
      BB3 = BB3 + t1
      
    }
    
  }
  
  -BB3
  
}
d2ldb32(rep(0,51), G2, k=16, cohort = 2, dummy)

d2ldag2 = function(theta, df, k = 16, g, cohort, dummy){ # validated
  # g controls groups
  
  b1 = theta[1:k]
  b2 = theta[(k+1):(2*k)]
  b3 = theta[(2*k+1):(3*k)]
  ag = theta[3*k+g]
  
  G = length(df)
  AAg = 0
  
  d = as.matrix(df[[g]])
  n = nrow(d)
  
  D = dummy[[cohort]]
  D1 = D[,1]; D2 = D[,2]; D3 = D[,3]; D4 = D[,4]
  
  if (G == 1){ range = 1:n } else {
    if (cohort == 1){ range = idx.d1[[g]] }
    if (cohort == 2){ range = idx.d2[[g]] }
  }
  D11 = D1[range]; D21 = D2[range]; D31 = D3[range]; D41 = D4[range]
  
  for (i in 1:n){
    x_i = as.vector(d[i,])
    eta1 = as.numeric(x_i %*% b1 + ag)
    eta2 = as.numeric(x_i %*% b2 + ag)
    eta3 = as.numeric(x_i %*% b3 + ag)
    
    t1 = plogis(eta1) * (1-plogis(eta1)) * D11[i]
    t2 = (plogis(eta1)*(1-plogis(eta1)) + plogis(eta2)*(1-plogis(eta2))) * D21[i]
    t3 = ((plogis(eta1)*(1-plogis(eta1)) + plogis(eta2)*(1-plogis(eta2)) + plogis(eta3)*(1-plogis(eta3)))) * (D31[i] + D41[i])
    
    AAg = AAg + (t1 + t2 + t3)
  }
  -AAg
}
d2ldag2(rep(0,51), df=G2, k=16, g=1, cohort=2, dummy)

d2ldb1dag = function(theta, df, k = 16, g, cohort, dummy){ # validated
  # g controls groups
  
  b1 = theta[1:k]
  b2 = theta[(k+1):(2*k)]
  b3 = theta[(2*k+1):(3*k)]
  ag = theta[3*k+g]
  
  G = length(df)
  B1A = rep(0, k)
  
  d = as.matrix(df[[g]])
  n = nrow(d)
  
  D = dummy[[cohort]]
  D1 = D[,1]; D2 = D[,2]; D3 = D[,3]; D4 = D[,4]
  
  if (G == 1){ range = 1:n } else {
    if (cohort == 1){ range = idx.d1[[g]] }
    if (cohort == 2){ range = idx.d2[[g]] }
  }
  D11 = D1[range]; D21 = D2[range]; D31 = D3[range]; D41 = D4[range]
  
  for (i in 1:n){
    x_i = as.vector(d[i,])
    eta1 = as.numeric(x_i %*% b1 + ag)
    #eta2 = as.numeric(x_i %*% b2 + ag)
    #eta3 = as.numeric(x_i %*% b3 + ag)
    
    t1 = plogis(eta1) * (1-plogis(eta1)) * x_i
    
    B1A = B1A + t1
  }
  -as.matrix(B1A, nr=k, nc=1)
}
d2ldb1dag(rep(0,50), G1, k = 16, g=1, cohort=1, dummy)
  
d2ldb2dag = function(theta, df, k = 16, g, cohort, dummy){ # validated
  # g controls groups
  
  b1 = theta[1:k]
  b2 = theta[(k+1):(2*k)]
  b3 = theta[(2*k+1):(3*k)]
  ag = theta[3*k+g]
  
  G = length(df)
  B2A = rep(0, k)
  d = as.matrix(df[[g]])
  n = nrow(d)
  
  D = dummy[[cohort]]
  D1 = D[,1]; D2 = D[,2]; D3 = D[,3]; D4 = D[,4]
  
  if (G == 1){ range = 1:n } else {
    if (cohort == 1){ range = idx.d1[[g]] }
    if (cohort == 2){ range = idx.d2[[g]] }
  }
  D11 = D1[range]; D21 = D2[range]; D31 = D3[range]; D41 = D4[range]
  
  for (i in 1:n){
    x_i = as.vector(d[i,])
    #eta1 = as.numeric(x_i %*% b1 + ag)
    eta2 = as.numeric(x_i %*% b2 + ag)
    #eta3 = as.numeric(x_i %*% b3 + ag)
    
    t1 = plogis(eta2) * (1-plogis(eta2)) * x_i * (D21[i] + D31[i] + D41[i])
    
    B2A = B2A + t1
  }
  -as.matrix(B2A, nr=k, nc=1)
}
d2ldb2dag(rep(0,50), G1, k = 16, g=2, cohort=1, dummy)

d2ldb3dag = function(theta, df, k = 16, g, cohort, dummy){ # validated
  # g controls groups
  
  b1 = theta[1:k]
  b2 = theta[(k+1):(2*k)]
  b3 = theta[(2*k+1):(3*k)]
  ag = theta[3*k+g]
  
  G = length(df)
  B3A = rep(0, k)
  d = as.matrix(df[[g]])
  n = nrow(d)
  
  D = dummy[[cohort]]
  D1 = D[,1]; D2 = D[,2]; D3 = D[,3]; D4 = D[,4]
  
  if (G == 1){ range = 1:n } else {
    if (cohort == 1){ range = idx.d1[[g]] }
    if (cohort == 2){ range = idx.d2[[g]] }
  }
  D11 = D1[range]; D21 = D2[range]; D31 = D3[range]; D41 = D4[range]
  
  for (i in 1:n){
    x_i = as.vector(d[i,])
    #eta1 = as.numeric(x_i %*% b1 + ag)
    #eta2 = as.numeric(x_i %*% b2 + ag)
    eta3 = as.numeric(x_i %*% b3 + ag)
    
    t1 = plogis(eta3) * (1-plogis(eta3)) * x_i * (D31[i] + D41[i])
    
    B3A = B3A + t1
  }
  -as.matrix(B3A, nr=k, nc=1)
}
d2ldb3dag(rep(0,50), G1, k = 16, g=1, cohort = 1, dummy)

d2ldagdb1 = function(theta, df, k = 16, g, cohort, dummy){ t(d2ldb1dag(theta, df, k = 16, g, cohort, dummy)) } # validated
d2ldagdb2 = function(theta, df, k = 16, g, cohort, dummy){ t(d2ldb2dag(theta, df, k = 16, g, cohort, dummy)) } # validated
d2ldagdb3 = function(theta, df, k = 16, g, cohort, dummy){ t(d2ldb3dag(theta, df, k = 16, g, cohort, dummy)) } # validated

hess = function(theta, df, k = 16, cohort, dummy){ # validated
  
  G = length(df)
  H = matrix(0, length(theta), length(theta))
  H[1:k,1:k] = d2ldb12(theta, df, k, dummy)
  H[((k+1):(2*k)),((k+1):(2*k))] = d2ldb22(theta, df, k, cohort, dummy)
  H[((2*k+1):(3*k)),((2*k+1):(3*k))] = d2ldb32(theta, df, k, cohort, dummy)
  
  for (g in 1:G){
    H[(3*k)+g, (3*k)+g] = d2ldag2(theta, df, k, g, cohort, dummy)
    
    H[1:k, (3*k)+g] = d2ldb1dag(theta, df, k, g, cohort, dummy)
    H[(k+1):(2*k), (3*k)+g] = d2ldb2dag(theta, df, k, g, cohort, dummy)
    H[(2*k+1):(3*k), (3*k)+g] = d2ldb3dag(theta, df, k, g, cohort, dummy)
    
    H[(3*k)+g, 1:k] = d2ldagdb1(theta, df, k, g, cohort, dummy)
    H[(3*k)+g, (k+1):(2*k)] = d2ldagdb2(theta, df, k, g, cohort, dummy)
    H[(3*k)+g, (2*k+1):(3*k)] = d2ldagdb3(theta, df, k, g, cohort, dummy)
  }
  H
  
}
hess(rep(0,49), k=16, list(X1.c), cohort = 1, dummy)
### End: Hessian Matrix ###

##### Implementation #####
require(pracma)
brd = function(theta, df, k = 16, cohort, dummy, maxiter = 100, tol = .Machine$double.eps^(1/2)) {
  # input df shall be a well-prepared data set, in the type of a list
  # standardize whatever it's needed for ESTIMATION
  x0 = theta
  loglik = ll(x0, df, k, cohort, dummy)
  print(loglik)
  
  y0 = score(x0, df, k, cohort, dummy)
  A0 = hess(x0, df, k, cohort, dummy)
  B0 = ginv(A0)
  
  xnew <- x0 - B0 %*% y0
  ynew <- score(xnew, df, k, cohort, dummy)
  niter <- 1
  while (niter < maxiter) {
    
    loglik = c(loglik, ll(xnew, df, k, cohort, dummy))
    print(ll(xnew, df, k, cohort, dummy))
    
    s <- xnew - x0
    d <- ynew - y0
    if (norm(s, "F") < tol || norm(as.matrix(ynew), "F") < tol) 
      break
    
    print(c(norm(s, "F"), norm(as.matrix(ynew), "F")))
    B0 <- B0 + (s - B0 %*% d) %*% t(s) %*% B0/c(t(s) %*% B0 %*% d)
    x0 <- xnew
    xnew <- xnew - B0 %*% ynew
    y0 <- ynew
    ynew <- score(xnew, df, k, cohort, dummy)
    niter <- niter + 1
  }
  if (niter >= maxiter) 
    warning(paste("Not converged: Max number of iterations reached."))
  fnew <- sqrt(sum(ynew^2))
  return(list(zero = c(xnew), fnorm = fnew, niter = niter, loglik = loglik, hess = hess(xnew, df, k, cohort, dummy)))
}