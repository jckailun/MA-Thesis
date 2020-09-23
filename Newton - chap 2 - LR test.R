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

ll.lr = function(theta, cohort, dummy){ # validated
  
  a = theta; G = length(a)
  l = 0
  D = dummy[[cohort]]
  D11 = D[,1]; D21 = D[,2]; D31 = D[,3]; D41 = D[,4]
  
  if ((cohort == 1) & (G != 1)) {n = list(idx.g1.d1, idx.g2.d1)} 
  if ((cohort == 2) & (G != 1)) {n = idx.d2} 
  if ((cohort == 1) & (G == 1)) {n = list(1:n1)} 
  if ((cohort == 2) & (G == 1)) {n = list(1:n2)} 
  
  for (j in 1:G){
    
    aj = a[j]

    for (i in n[[j]]){
      
      eta1 = as.numeric(aj)
      eta2 = as.numeric(aj)
      eta3 = as.numeric(aj)
      
      t1 = -log(1 + exp(eta1))*D11[i]
      t2 = (eta1 - log(1 + exp(eta1)) - log(1 + exp(eta2))) * D21[i]
      t3 = (eta1 + eta2 - log(1 + exp(eta1)) - log(1 + exp(eta2)) - log(1 + exp(eta3))) * D31[i]
      t4 = (eta1 + eta2 + eta3 - log(1 + exp(eta1)) - log(1 + exp(eta2)) - log(1 + exp(eta3))) * D41[i]
      
      l = l + (t1 + t2 + t3 + t4)
      
    }
  }
  l
}
ll.lr(0, cohort=1, dummy.bs)
### Score Function ###

dldag.lr = function(theta, g, cohort, dummy){ # validated
  # g controls groups
  
  a = theta; G = length(a)
  l = 0
  D = dummy[[cohort]]
  D11 = D[,1]; D21 = D[,2]; D31 = D[,3]; D41 = D[,4]
  
  if ((cohort == 1) & (G != 1)) {n = list(idx.g1.d1, idx.g2.d1)} 
  if ((cohort == 2) & (G != 1)) {n = idx.d2} 
  if ((cohort == 1) & (G == 1)) {n = list(1:n1)} 
  if ((cohort == 2) & (G == 1)) {n = list(1:n2)} 
  
  Ag = 0
  aj = a[g]

  for (i in n[[g]]){
    
    eta1 = as.numeric(aj)
    eta2 = as.numeric(aj)
    eta3 = as.numeric(aj)
    
    t1 = - plogis(eta1) * D11[i]
    t2 = (1 - plogis(eta1) - plogis(eta2))  * D21[i]
    t3 = (2 - plogis(eta1) - plogis(eta2) - plogis(eta3)) * D31[i]
    t4 = (3 - plogis(eta1) - plogis(eta2) - plogis(eta3)) * D41[i]
    
    Ag = Ag + (t1 + t2 + t3 + t4)
  }
  Ag
}
dldag.lr(0, g=1, cohort=1, dummy.bs)

score.lr = function(theta, cohort, dummy){ # validated
  
  G = length(theta)
  sa = NULL
  for (j in 1:G){
    sa = c(sa, dldag.lr(theta, g = j, cohort, dummy))
  }
  
  sa
  
}
score.lr(0, cohort = 1, dummy.bs)
### End: Score Function ###
### Hessian Matrix ###

d2ldag2.lr = function(theta, g, cohort, dummy){ # validated
  # g controls groups
  ag = theta
  
  G = length(ag)
  D = dummy[[cohort]]
  D11 = D[,1]; D21 = D[,2]; D31 = D[,3]; D41 = D[,4]
  
  AAg = 0
  if ((cohort == 1) & (G != 1)) {n = list(idx.g1.d1, idx.g2.d1)} 
  if ((cohort == 2) & (G != 1)) {n = idx.d2} 
  if ((cohort == 1) & (G == 1)) {n = list(1:n1)} 
  if ((cohort == 2) & (G == 1)) {n = list(1:n2)} 

  for (i in n[[g]]){
    
    eta1 = as.numeric(ag)
    eta2 = as.numeric(ag)
    eta3 = as.numeric(ag)
    
    t1 = plogis(eta1) * (1-plogis(eta1)) * D11[i]
    t2 = (plogis(eta1) * (1-plogis(eta1)) + plogis(eta2)*(1-plogis(eta2))) * D21[i]
    t3 = ((plogis(eta1) * (1-plogis(eta1)) + plogis(eta2)*(1-plogis(eta2)) + plogis(eta3)*(1-plogis(eta3)))) * (D31[i] + D41[i])
    
    AAg = AAg + (t1 + t2 + t3)
  }
  -AAg
}
d2ldag2.lr(0, g=1, cohort=1, dummy.bs)

hess.lr = function(theta, cohort, dummy){ # validated
  
  G = length(theta)
  H = matrix(0, length(theta), length(theta))
  
  for (g in 1:G){
    H[g,g] = d2ldag2.lr(theta, g, cohort, dummy)[g]
  }
  H
}
hess.lr(c(0,0), cohort = 2, dummy.bs)
### End: Hessian Matrix ###

##### Implementation #####
require(pracma)
brd.lr = function(theta, cohort, dummy, maxiter = 100, tol = .Machine$double.eps^(1/2)) {
  # input df shall be a well-prepared data set
  # standardize whatever it's needed for ESTIMATION
  x0 = theta
  loglik = ll.lr(x0, cohort, dummy)
  print(loglik)
  
  y0 = score.lr(x0, cohort, dummy)
  A0 = hess.lr(x0, cohort, dummy)
  B0 = inv(A0)
  
  xnew <- x0 - B0 %*% y0
  ynew <- score.lr(xnew, cohort, dummy)
  niter <- 1
  while (niter < maxiter) {
    
    loglik = c(loglik, ll.lr(xnew, cohort, dummy))
    print(ll.lr(xnew, cohort, dummy))
    
    s <- xnew - x0
    d <- ynew - y0
    if (norm(s, "F") < tol || norm(as.matrix(ynew), "F") < tol) 
      break
    
    print(c(norm(s, "F"), norm(as.matrix(ynew), "F")))
    B0 <- B0 + (s - B0 %*% d) %*% t(s) %*% B0/c(t(s) %*% B0 %*% d)
    x0 <- xnew
    xnew <- xnew - B0 %*% ynew
    y0 <- ynew
    ynew <- score.lr(xnew, cohort, dummy)
    niter <- niter + 1
  }
  if (niter >= maxiter) 
    warning(paste("Not converged: Max number of iterations reached."))
  fnew <- sqrt(sum(ynew^2))
  return(list(zero = c(xnew), fnorm = fnew, niter = niter, loglik = loglik, hess = hess.lr(xnew, cohort, dummy)))
}