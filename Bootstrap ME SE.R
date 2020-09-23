### Bootstrap ME s.e. ###
### Cohort 1 ###
B = 1000
set.seed(1000)
est01.bs = NULL; ME.bs = list(matrix(0, nc=16, nr=B), matrix(0, nc=16, nr=B), matrix(0, nc=16, nr=B))
for (b in 1:B){
  
  idx.bs = sample(1:n1, n1, replace = TRUE)
  dummy.bs = dummy
  dummy.bs[[1]] = dummy.bs[[1]][idx.bs, ]
  
  theta.01 = rep(0, 49)
  X1.bs = X1.c[idx.bs, ]
  out = brd(theta.01, list(X1.bs), k=16, cohort = 1, dummy.bs)
  est01.bs = rbind(est01.bs, out$zero)
  
  b1.bs = list(out$zero[1:16], out$zero[17:32], out$zero[33:48])
  a.bs = out$zero[49]
  
  for (j in 1:3){
    for (k in 1:16){
      if (k %in% c(11, 13:16)){
        
        dat.temp = X1.bs
        
        dat.temp[, k] = 1
        eta1.bs = dat.temp %*% b1.bs[[j]] + a.bs
        
        dat.temp[, k] = 0
        eta0.bs = dat.temp %*% b1.bs[[j]] + a.bs
        
        ME.bs[[j]][b, k] = mean(plogis(eta1.bs) - plogis(eta0.bs))
        
      } else {
        
        eta.bs = X1.bs %*% b1.bs[[j]] + a.bs
        ME.bs[[j]][b, k] = mean( plogis(eta.bs) * (1-plogis(eta.bs)) * b1.bs[[j]][k])
    
      }
    }
  }
}
se = NULL
for (j in 1:3){
  se = cbind(se, sqrt(diag(cov(ME.bs[[j]]))))
}
    
ME1.out = cbind(J[,1], se[,1], J[,1]/se[,1],
               J[,2], se[,2], J[,2]/se[,2],
               J[,3], se[,3], J[,3]/se[,3])
write.csv(ME1.out, 'Marginal Effects of Cohort 1.csv')
    

### Cohort 2 ###

B = 1000
set.seed(1000)
est02.bs = NULL; ME2.bs = list(matrix(0, nc=16, nr=B), matrix(0, nc=16, nr=B), matrix(0, nc=16, nr=B))
for (b in 1:B){
  
  idx.bs = sample(1:n2, n2, replace = TRUE)
  dummy.bs = dummy
  dummy.bs[[2]] = dummy.bs[[2]][idx.bs, ]
  
  theta.02 = rep(0, 49)
  X2.bs = X2.c[idx.bs, ]
  out = brd(theta.02, list(X2.bs), k=16, cohort = 2, dummy.bs)
  est02.bs = rbind(est02.bs, out$zero)
  
  b2.bs = list(out$zero[1:16], out$zero[17:32], out$zero[33:48])
  a2.bs = out$zero[49]
  
  for (j in 1:3){
    for (k in 1:16){
      if (k %in% c(11, 13:16)){
        
        dat.temp = X2.bs
        
        dat.temp[, k] = 1
        eta1.bs = dat.temp %*% b2.bs[[j]] + a2.bs
        
        dat.temp[, k] = 0
        eta0.bs = dat.temp %*% b2.bs[[j]] + a2.bs
        
        ME2.bs[[j]][b, k] = mean(plogis(eta1.bs) - plogis(eta0.bs))
        
      } else {
        
        eta.bs = X2.bs %*% b2.bs[[j]] + a2.bs
        ME2.bs[[j]][b, k] = mean( plogis(eta.bs) * (1-plogis(eta.bs)) * b2.bs[[j]][k])
        
      }
    }
  }
}
se2 = NULL
for (j in 1:3){
  se2 = cbind(se2, sqrt(diag(cov(ME2.bs[[j]]))))
}

ME2.out = cbind(J2[,1], se2[,1], J2[,1]/se2[,1],
                J2[,2], se2[,2], J2[,2]/se2[,2],
                J2[,3], se2[,3], J2[,3]/se2[,3])
write.csv(ME2.out, 'Marginal Effects of Cohort 2.csv')

