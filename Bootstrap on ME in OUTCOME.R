B = 1000
set.seed(1000)
est01.bs = NULL
ME.bs.c = list(matrix(0, nc=16, nr=B), matrix(0, nc=16, nr=B), matrix(0, nc=16, nr=B), matrix(0, nc=16, nr=B))
for (b in 1:B){
  
  idx.bs = sample(1:n1, n1, replace = TRUE)
  dummy.bs = dummy
  dummy.bs[[1]] = dummy.bs[[1]][idx.bs, ]
  
  theta.01 = rep(0, 49)
  X1.bs = X1.c[idx.bs, ]
  out = brd(theta.01, list(X1.bs), k=16, cohort = 1, dummy.bs)
  est01.bs = rbind(est01.bs, out$zero)
  
  b1.bs = list(out$zero[1:16], out$zero[17:32], out$zero[33:48])
  a1.bs = out$zero[49]
  
    for (k in 1:16){
      if (k %in% c(11, 13:16)){
        
        dat.temp = X1.bs
        
        dat.temp[, k] = 1
        eta1.bs = dat.temp %*% b1.bs[[1]] + a1.bs
        eta12.bs = dat.temp %*% b1.bs[[2]]+ a1.bs
        eta13.bs = dat.temp %*% b1.bs[[3]]+ a1.bs
        
        dat.temp[, k] = 0
        eta0.bs = dat.temp %*% b1.bs[[1]]+ a1.bs
        eta02.bs = dat.temp %*% b1.bs[[2]]+ a1.bs
        eta03.bs = dat.temp %*% b1.bs[[3]]+ a1.bs
        
        #ME.bs.c[1][b, k] = mean(plogis(eta1.bs) - plogis(eta0.bs))
        
        ME.bs.c[[1]][b, k] = mean((1-plogis(eta1.bs)) - (1-plogis(eta0.bs)))
        
        ME.bs.c[[2]][b, k] = mean(plogis(eta1.bs)*(1-plogis(eta12.bs)) - plogis(eta0.bs)*(1-plogis(eta02.bs)))
        
        ME.bs.c[[3]][b, k] = mean(plogis(eta1.bs)*(plogis(eta12.bs)*(1-plogis(eta13.bs))) - plogis(eta0.bs)*plogis(eta02.bs)*(1-plogis(eta03.bs)))
        
        ME.bs.c[[4]][b, k] = mean(plogis(eta1.bs)*(plogis(eta12.bs)*plogis(eta13.bs)) - plogis(eta0.bs)*plogis(eta02.bs)*plogis(eta03.bs))
        
      } else {
        
        eta1.bs = X1.bs %*% b1.bs[[1]]+ a1.bs
        eta2.bs = X1.bs %*% b1.bs[[2]]+ a1.bs
        eta3.bs = X1.bs %*% b1.bs[[3]]+ a1.bs
        
        ME.bs.c[[1]][b, k] = mean(  - plogis(eta1.bs) * (1-plogis(eta1.bs)) * b1.bs[[1]][k] )
        ME.bs.c[[2]][b, k] = mean(  (1 - plogis(eta2.bs)) * plogis(eta1.bs) * (1-plogis(eta1.bs)) *  b1.bs[[1]][k] +
                                   plogis(eta1.bs) * (- plogis(eta2.bs) * (1-plogis(eta2.bs))) * b1.bs[[2]][k]  ) 
        ME.bs.c[[3]][b, k] = mean( (1 - plogis(eta3.bs)) * plogis(eta2.bs) * (plogis(eta1.bs) * (1-plogis(eta1.bs))) * b1.bs[[1]][k] +
                                  (1 - plogis(eta3.bs)) * plogis(eta1.bs) * (plogis(eta2.bs) * (1-plogis(eta2.bs))) * b1.bs[[2]][k] + 
                                  plogis(eta1.bs) * plogis(eta2.bs) * (-plogis(eta3.bs) * (1-plogis(eta3.bs))) * b1.bs[[3]][k] )
        ME.bs.c[[4]][b, k] = mean( plogis(eta3.bs) * plogis(eta2.bs) * plogis(eta1.bs) * (1-plogis(eta1.bs)) * b1.bs[[1]][k] +
                                   plogis(eta3.bs) * plogis(eta1.bs) * plogis(eta2.bs) * (1-plogis(eta2.bs)) * b1.bs[[2]][k] + 
                                   plogis(eta1.bs) * plogis(eta2.bs) * plogis(eta3.bs) * (1-plogis(eta3.bs)) * b1.bs[[3]][k] ) 
        
      }
    }}

se.c1 = NULL
for (j in 1:4){
  se.c1 = cbind(se.c1, sqrt(diag(cov(ME.bs.c[[j]]))))
}

ME1.class[,2] = se.c1[,1]; ME1.class[,3] = ME1.class[,1] / ME1.class[,2]
ME1.class[,5] = se.c1[,2]; ME1.class[,6] = ME1.class[,4] / ME1.class[,5]
ME1.class[,8] = se.c1[,3]; ME1.class[,9] = ME1.class[,7] / ME1.class[,8]
ME1.class[,11] = se.c1[,4]; ME1.class[,12] = ME1.class[,10] / ME1.class[,11]

pval.me1 = matrix(0, nr=16, nc=4)
ref = ME1.class[,c(1,4,7,10)]
for (j in 1:4){
  for (k in 1:16){
    
    pval.me1[k, j] = mean(abs(ME.bs.c[[j]][,k]-ref[,j][k])>=abs(ref[,j][k]))
  }}


write.csv(cbind(ME1.class, pval.me1), "Marginal Effects on OUTCOME - Cohort 1.csv")
### Cohort 2 ###

B = 1000
set.seed(1000)
est01.bs = NULL
ME2.bs.c = list(matrix(0, nc=16, nr=B), matrix(0, nc=16, nr=B), matrix(0, nc=16, nr=B), matrix(0, nc=16, nr=B))
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
  
  for (k in 1:16){
    if (k %in% c(11, 13:16)){
      
      dat.temp = X2.bs
      
      dat.temp[, k] = 1
      eta1.bs = dat.temp %*% b2.bs[[1]] + a2.bs
      eta12.bs = dat.temp %*% b2.bs[[2]]+ a2.bs
      eta13.bs = dat.temp %*% b2.bs[[3]]+ a2.bs
      
      dat.temp[, k] = 0
      eta0.bs = dat.temp %*% b2.bs[[1]]+ a2.bs
      eta02.bs = dat.temp %*% b2.bs[[2]]+ a2.bs
      eta03.bs = dat.temp %*% b2.bs[[3]]+ a2.bs
      
      #ME.bs.c[1][b, k] = mean(plogis(eta1.bs) - plogis(eta0.bs))
      
      ME2.bs.c[[1]][b, k] = mean((1-plogis(eta1.bs)) - (1-plogis(eta0.bs)))
      
      ME2.bs.c[[2]][b, k] = mean(plogis(eta1.bs)*(1-plogis(eta12.bs)) - plogis(eta0.bs)*(1-plogis(eta02.bs)))
      
      ME2.bs.c[[3]][b, k] = mean(plogis(eta1.bs)*(plogis(eta12.bs)*(1-plogis(eta13.bs))) - plogis(eta0.bs)*plogis(eta02.bs)*(1-plogis(eta03.bs)))
      
      ME2.bs.c[[4]][b, k] = mean(plogis(eta1.bs)*(plogis(eta12.bs)*plogis(eta13.bs)) - plogis(eta0.bs)*plogis(eta02.bs)*plogis(eta03.bs))
      
    } else {
      
      eta1.bs = X2.bs %*% b2.bs[[1]]+ a2.bs
      eta2.bs = X2.bs %*% b2.bs[[2]]+ a2.bs
      eta3.bs = X2.bs %*% b2.bs[[3]]+ a2.bs
      
      ME2.bs.c[[1]][b, k] = mean(  - plogis(eta1.bs) * (1-plogis(eta1.bs)) * b2.bs[[1]][k] )
      ME2.bs.c[[2]][b, k] = mean(  (1 - plogis(eta2.bs)) * plogis(eta1.bs) * (1-plogis(eta1.bs)) *  b2.bs[[1]][k] +
                                    plogis(eta1.bs) * (- plogis(eta2.bs) * (1-plogis(eta2.bs))) * b2.bs[[2]][k]  ) 
      ME2.bs.c[[3]][b, k] = mean( (1 - plogis(eta3.bs)) * plogis(eta2.bs) * (plogis(eta1.bs) * (1-plogis(eta1.bs))) * b2.bs[[1]][k] +
                                   (1 - plogis(eta3.bs)) * plogis(eta1.bs) * (plogis(eta2.bs) * (1-plogis(eta2.bs))) * b2.bs[[2]][k] + 
                                   plogis(eta1.bs) * plogis(eta2.bs) * (-plogis(eta3.bs) * (1-plogis(eta3.bs))) * b2.bs[[3]][k] )
      ME2.bs.c[[4]][b, k] = mean( plogis(eta3.bs) * plogis(eta2.bs) * plogis(eta1.bs) * (1-plogis(eta1.bs)) * b2.bs[[1]][k] +
                                   plogis(eta3.bs) * plogis(eta1.bs) * plogis(eta2.bs) * (1-plogis(eta2.bs)) * b2.bs[[2]][k] + 
                                   plogis(eta1.bs) * plogis(eta2.bs) * plogis(eta3.bs) * (1-plogis(eta3.bs)) * b2.bs[[3]][k] ) 
      
    }
  }}

se.c2 = NULL
for (j in 1:4){
  se.c2 = cbind(se.c2, sqrt(diag(cov(ME2.bs.c[[j]]))))
}

ME2.class[,2] = se.c2[,1]; ME2.class[,3] = ME2.class[,1] / ME2.class[,2]
ME2.class[,5] = se.c2[,2]; ME2.class[,6] = ME2.class[,4] / ME2.class[,5]
ME2.class[,8] = se.c2[,3]; ME2.class[,9] = ME2.class[,7] / ME2.class[,8]
ME2.class[,11] = se.c2[,4]; ME2.class[,12] = ME2.class[,10] / ME2.class[,11]


pval.me2 = matrix(0, nr=16, nc=4)
ref = ME2.class[,c(1,4,7,10)]
for (j in 1:4){
for (k in 1:16){
 
pval.me2[k, j] = mean(abs(ME2.bs.c[[j]][,k]-ref[,j][k])>=abs(ref[,j][k]))
}}
write.csv(cbind(ME2.class, pval.me2), "Marginal Effects on OUTCOME - Cohort 2.csv")

