##### Newton's method for CV on Cohort 1 #####

ll.cv = function(theta, df, k = 16, cohort, id){ # validated
  
  G = length(df)
  b1 = theta[1:k]
  b2 = theta[(k+1):(2*k)]
  b3 = theta[(2*k+1):(3*k)]
  a = theta[-(1:(3*k))]
  l = 0
  
  for (j in 1:G){
    
    d = as.matrix(df[[j]])
    aj = a[j]
    n = nrow(d)
    
    D = dummy[[cohort]][id[[j]],]
    D11 = D[,1]; D21 = D[,2]; D31 = D[,3]; D41 = D[,4]
    
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
ll.cv(rep(0, 50), train.G1, cohort = 1, k=16, id = index.cv.11)
### Score Function ###
dldb1.cv = function(theta, df, k = 16, cohort, id){ # validated
  
  G = length(df)
  b1 = theta[1:k]
  b2 = theta[(k+1):(2*k)]
  b3 = theta[(2*k+1):(3*k)]
  a = theta[-(1:(3*k))]
  B1 = 0
  
  for (j in 1:G){
    
    d = as.matrix(df[[j]])
    aj = a[j]
    n = nrow(d)
    
    D = dummy[[cohort]][id[[j]],]
    D11 = D[,1]; D21 = D[,2]; D31 = D[,3]; D41 = D[,4]
    
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
dldb1.cv(rep(0,50), df=train.G1, k=16, cohort=1, id = index.cv.11)

dldb2.cv = function(theta, df, k=16, cohort, id){ 
  
  G = length(df)
  b1 = theta[1:k]
  b2 = theta[(k+1):(2*k)]
  b3 = theta[(2*k+1):(3*k)]
  a = theta[-(1:(3*k))]
  B2 = 0
  
  for (j in 1:G){
    
    d = as.matrix(df[[j]])
    aj = a[j]
    n = nrow(d)
    
    D = dummy[[cohort]][id[[j]],]
    D11 = D[,1]; D21 = D[,2]; D31 = D[,3]; D41 = D[,4]
    
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
dldb2.cv(rep(0,50), df=train.G1, k=16, cohort=1, id =index.cv.11)

dldb3.cv = function(theta, df, k=16, cohort, id){ # validated
  
  G = length(df)
  b1 = theta[1:k]
  b2 = theta[(k+1):(2*k)]
  b3 = theta[(2*k+1):(3*k)]
  a = theta[-(1:(3*k))]
  B3 = 0
  
  for (j in 1:G){
    
    d = as.matrix(df[[j]])
    aj = a[j]
    n = nrow(d)
    
    D = dummy[[cohort]][id[[j]],]
    D11 = D[,1]; D21 = D[,2]; D31 = D[,3]; D41 = D[,4]
    
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
dldb3.cv(rep(0,50), df=train.G1, k=16, cohort=1, id = index.cv.11)

dldag.cv = function(theta, df, k = 16, g, cohort, id){ # validated
  # g controls groups
  
  b1 = theta[1:k]
  b2 = theta[(k+1):(2*k)]
  b3 = theta[(2*k+1):(3*k)]
  aj = theta[3*k+g]
  
  G = length(df)
  D = dummy[[cohort]][id[[g]],]
  D11 = D[,1]; D21 = D[,2]; D31 = D[,3]; D41 = D[,4]
  
  Ag = 0
  
  d = as.matrix(df[[g]])
  n = nrow(d)
  
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
dldag.cv(rep(0,50), df=train.G1, k=16, g=1, cohort=1, id = index.cv.11)

score.cv = function(theta, df, k=16, cohort, id){ # validated
  
  G = length(df)
  sa = NULL
  
  sb1 = dldb1.cv(theta, df, k = 16, cohort, id)
  sb2 = dldb2.cv(theta, df, k = 16, cohort, id)
  sb3 = dldb3.cv(theta, df, k = 16, cohort, id)
  
  for (j in 1:G){
    sa = c(sa, dldag.cv(theta, df, k = 16, g = j, cohort, id))
  }
  
  c(sb1, sb2, sb3, sa)
  
}
score.cv(rep(0,50), train.G1, cohort=1, id=index.cv.11)
### End: Score Function ###

### Hessian Matrix ###
d2ldb12.cv = function(theta, df, k = 16){ # validated
  
  G = length(df)
  b1 = theta[1:k]
  b2 = theta[(k+1):(2*k)]
  b3 = theta[(2*k+1):(3*k)]
  a = theta[-(1:(3*k))]
  
  BB1 = matrix(0, nr = k, nc = k)
  
  for (g in 1:G){
    
    d = as.matrix(df[[g]])
    ag = a[g]
    n = nrow(d)
    
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
d2ldb12.cv(rep(0,50), train.G1)

d2ldb22.cv = function(theta, df, k = 16, cohort, id){ # validated
  
  G = length(df)
  b1 = theta[1:k]
  b2 = theta[(k+1):(2*k)]
  b3 = theta[(2*k+1):(3*k)]
  a = theta[-(1:(3*k))]
  
  BB2 = matrix(0, nr = k, nc = k)
  
  for (g in 1:G){
    
    d = as.matrix(df[[g]])
    ag = a[g]
    n = nrow(d)
    D = dummy[[cohort]][id[[g]],]
    D11 = D[,1]; D21 = D[,2]; D31 = D[,3]; D41 = D[,4]
    
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
d2ldb22.cv(rep(0,50), train.G1, cohort = 1, id=index.cv.11)

d2ldb32.cv = function(theta, df, k = 16, cohort, id){ # validated
  
  G = length(df)
  b1 = theta[1:k]
  b2 = theta[(k+1):(2*k)]
  b3 = theta[(2*k+1):(3*k)]
  a = theta[-(1:(3*k))]
  
  BB3 = matrix(0, nr = k, nc = k)
 
  for (g in 1:G){
    
    d = as.matrix(df[[g]])
    ag = a[g]
    n = nrow(d)
    D = dummy[[cohort]][id[[g]],]
    D11 = D[,1]; D21 = D[,2]; D31 = D[,3]; D41 = D[,4]
    
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
d2ldb32.cv(rep(0,50), train.G1, cohort = 1, id=index.cv.11)

d2ldag2.cv = function(theta, df, k = 16, g, cohort, id){ # validated
  # g controls groups
  
  b1 = theta[1:k]
  b2 = theta[(k+1):(2*k)]
  b3 = theta[(2*k+1):(3*k)]
  ag = theta[3*k+g]
  
  G = length(df)
  D = dummy[[cohort]][id[[g]],]
  D11 = D[,1]; D21 = D[,2]; D31 = D[,3]; D41 = D[,4]
  
  AAg = 0
  
  d = as.matrix(df[[g]])
  n = nrow(d)
  
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
d2ldag2.cv(rep(0,50), df=train.G1, k=16, g=1, cohort=1, index.cv.11)

d2ldb1dag.cv = function(theta, df, k = 16, g, cohort, id){ # validated
  # g controls groups
  
  b1 = theta[1:k]
  b2 = theta[(k+1):(2*k)]
  b3 = theta[(2*k+1):(3*k)]
  ag = theta[3*k+g]
  
  G = length(df)
  D = dummy[[cohort]][id[[g]],]
  D11 = D[,1]; D21 = D[,2]; D31 = D[,3]; D41 = D[,4]
  B1A = rep(0, k)
  
  d = as.matrix(df[[g]])
  n = nrow(d)
  
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
d2ldb1dag.cv(rep(0,50), train.G1, k = 16, g=1, cohort=1, index.cv.11)

d2ldb2dag.cv = function(theta, df, k = 16, g, cohort, id){ # validated
  # g controls groups
  
  b1 = theta[1:k]
  b2 = theta[(k+1):(2*k)]
  b3 = theta[(2*k+1):(3*k)]
  ag = theta[3*k+g]
  
  G = length(df)
  D = dummy[[cohort]][id[[g]],]
  D11 = D[,1]; D21 = D[,2]; D31 = D[,3]; D41 = D[,4]
  
  B2A = rep(0, k)
  
  d = as.matrix(df[[g]])
  n = nrow(d)
  
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
d2ldb2dag.cv(rep(0,50), train.G1, k = 16, g=1, cohort=1, index.cv.11)

d2ldb3dag.cv = function(theta, df, k = 16, g, cohort, id){ # validated
  # g controls groups
  
  b1 = theta[1:k]
  b2 = theta[(k+1):(2*k)]
  b3 = theta[(2*k+1):(3*k)]
  ag = theta[3*k+g]
  
  G = length(df)
  D = dummy[[cohort]][id[[g]],]
  D11 = D[,1]; D21 = D[,2]; D31 = D[,3]; D41 = D[,4]
  
  B3A = rep(0, k)
  
  d = as.matrix(df[[g]])
  n = nrow(d)
  
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
d2ldb3dag.cv(rep(0,50), train.G1, k = 16, g=1, cohort = 2, index.cv.11)

d2ldagdb1.cv = function(theta, df, k = 16, g, cohort, id){ t(d2ldb1dag.cv(theta, df, k = 16, g, cohort, id)) } # validated
d2ldagdb2.cv = function(theta, df, k = 16, g, cohort, id){ t(d2ldb2dag.cv(theta, df, k = 16, g, cohort, id)) } # validated
d2ldagdb3.cv = function(theta, df, k = 16, g, cohort, id){ t(d2ldb3dag.cv(theta, df, k = 16, g, cohort, id)) } # validated

hess.cv = function(theta, df, k = 16, cohort, id){ # validated
  
  G = length(df)
  H = matrix(0, length(theta), length(theta))
  H[1:k,1:k] = d2ldb12.cv(theta, df, k)
  H[((k+1):(2*k)),((k+1):(2*k))] = d2ldb22.cv(theta, df, k, cohort, id)
  H[((2*k+1):(3*k)),((2*k+1):(3*k))] = d2ldb32.cv(theta, df, k, cohort, id)
  
  for (g in 1:G){
    H[(3*k)+g, (3*k)+g] = d2ldag2.cv(theta, df, k, g, cohort, id)
    
    H[1:k, (3*k)+g] = d2ldb1dag.cv(theta, df, k, g, cohort, id)
    H[(k+1):(2*k), (3*k)+g] = d2ldb2dag.cv(theta, df, k, g, cohort, id)
    H[(2*k+1):(3*k), (3*k)+g] = d2ldb3dag.cv(theta, df, k, g, cohort, id)
    
    H[(3*k)+g, 1:k] = d2ldagdb1.cv(theta, df, k, g, cohort, id)
    H[(3*k)+g, (k+1):(2*k)] = d2ldagdb2.cv(theta, df, k, g, cohort, id)
    H[(3*k)+g, (2*k+1):(3*k)] = d2ldagdb3.cv(theta, df, k, g, cohort, id)
  }
  H
  
}
hess.cv(rep(0,50), train.G1, cohort = 1, id = index.cv.11)
### End: Hessian Matrix ###

##### Implementation #####
require(pracma)
brd.cv = function(theta, df, k = 16, cohort, id, maxiter = 100, tol = .Machine$double.eps^(1/2)) {
  # input df shall be a well-prepared data set
  # standardize whatever it's needed for ESTIMATION
  x0 = theta
  loglik = ll.cv(x0, df, k, cohort, id)
  print(loglik)
  
  y0 = score.cv(x0, df, k, cohort, id)
  A0 = hess.cv(x0, df, k, cohort, id)
  B0 = inv(A0)
  
  xnew <- x0 - B0 %*% y0
  ynew <- score.cv(xnew, df, k, cohort, id)
  niter <- 1
  while (niter < maxiter) {
    
    loglik = c(loglik, ll.cv(xnew, df, k, cohort, id))
    print(ll.cv(xnew, df, k, cohort, id))
    
    s <- xnew - x0
    d <- ynew - y0
    if (norm(s, "F") < tol || norm(as.matrix(ynew), "F") < tol) 
      break
    
    print(c(norm(s, "F"), norm(as.matrix(ynew), "F")))
    B0 <- B0 + (s - B0 %*% d) %*% t(s) %*% B0/c(t(s) %*% B0 %*% d)
    x0 <- xnew
    xnew <- xnew - B0 %*% ynew
    y0 <- ynew
    ynew <- score.cv(xnew, df, k, cohort, id)
    niter <- niter + 1
  }
  if (niter >= maxiter) 
    warning(paste("Not converged: Max number of iterations reached."))
  fnew <- sqrt(sum(ynew^2))
  return(list(zero = c(xnew), fnorm = fnew, niter = niter, loglik = loglik, hess = hess.cv(xnew, df, k, cohort, id)))
}