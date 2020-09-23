##### Marginal Effects #####
b1 = est01$zero; a1 = b1[49]
b1 = list(b1[1:16], b1[17:32], b1[33:48]) 
ME1.class = matrix(0, nc = 12, nr = 16)
idj = list(1:3, 4:6, 7:9, 10:12)

for (k in 1:16){
  
  if (k %in% c(11, 13:16)){ 

    dat.temp = X1.c
      
    dat.temp[, k] = 1
    eta1 = dat.temp %*% b1[[1]] + a1
    eta12 = dat.temp %*% b1[[2]]+ a1
    eta13 = dat.temp %*% b1[[3]]+ a1
      
    dat.temp[, k] = 0
    eta0 = dat.temp %*% b1[[1]]+ a1
    eta02 = dat.temp %*% b1[[2]]+ a1
    eta03 = dat.temp %*% b1[[3]]+ a1

    ME1.class[k, 1] = mean((1-plogis(eta1)) - (1-plogis(eta0)))

    ME1.class[k, 4] = mean(plogis(eta1)*(1-plogis(eta12)) - plogis(eta0)*(1-plogis(eta02)))

    ME1.class[k, 7] = mean(plogis(eta1)*(plogis(eta12)*(1-plogis(eta13))) - plogis(eta0)*plogis(eta02)*(1-plogis(eta03)))

    ME1.class[k, 10] = mean(plogis(eta1)*(plogis(eta12)*plogis(eta13)) - plogis(eta0)*plogis(eta02)*plogis(eta03))

    } else {
      
    eta1 = X1.c %*% b1[[1]]+ a1
    eta2 = X1.c %*% b1[[2]]+ a1
    eta3 = X1.c %*% b1[[3]]+ a1
    
    ME1.class[k, 1] = mean(  - plogis(eta1) * (1-plogis(eta1)) * b1[[1]][k] )
    ME1.class[k, 4] = mean(  (1 - plogis(eta2)) * plogis(eta1) * (1-plogis(eta1)) *  b1[[1]][k] +
                         plogis(eta1) * (- plogis(eta2) * (1-plogis(eta2))) * b1[[2]][k]  ) 
    ME1.class[k, 7] = mean( (1 - plogis(eta3)) * plogis(eta2) * (plogis(eta1) * (1-plogis(eta1))) * b1[[1]][k] +
                          (1 - plogis(eta3)) * plogis(eta1) * (plogis(eta2) * (1-plogis(eta2))) * b1[[2]][k] + 
                          plogis(eta1) * plogis(eta2) * (-plogis(eta3) * (1-plogis(eta3))) * b1[[3]][k] )
    ME1.class[k, 10] = mean( plogis(eta3) * plogis(eta2) * plogis(eta1) * (1-plogis(eta1)) * b1[[1]][k] +
                          plogis(eta3) * plogis(eta1) * plogis(eta2) * (1-plogis(eta2)) * b1[[2]][k] + 
                          plogis(eta1) * plogis(eta2) * plogis(eta3) * (1-plogis(eta3)) * b1[[3]][k] ) 
 
    }
}
 


### Cohort 2 ###
b2 = est02$zero; a2 = b2[49]
b2 = list(b2[1:16], b2[17:32], b2[33:48])
ME2.class = matrix(0, nc = 12, nr = 16)


for (k in 1:16){
  
  if (k %in% c(11, 13:16)){ 
    
    dat.temp = X2.c
    
    dat.temp[, k] = 1
    eta1 = dat.temp %*% b2[[1]] + a2
    eta12 = dat.temp %*% b2[[2]]+ a2
    eta13 = dat.temp %*% b2[[3]]+ a2
    
    dat.temp[, k] = 0
    eta0 = dat.temp %*% b2[[1]]+ a2
    eta02 = dat.temp %*% b2[[2]]+ a2
    eta03 = dat.temp %*% b2[[3]]+ a2
    
    ME2.class[k, 1] = mean((1-plogis(eta1)) - (1-plogis(eta0)))
    
    ME2.class[k, 4] = mean(plogis(eta1)*(1-plogis(eta12)) - plogis(eta0)*(1-plogis(eta02)))
    
    ME2.class[k, 7] = mean(plogis(eta1)*(plogis(eta12)*(1-plogis(eta13))) - plogis(eta0)*plogis(eta02)*(1-plogis(eta03)))
    
    ME2.class[k, 10] = mean(plogis(eta1)*(plogis(eta12)*plogis(eta13)) - plogis(eta0)*plogis(eta02)*plogis(eta03))
    
  } else {
    
    eta1 = X2.c %*% b2[[1]]+ a2
    eta2 = X2.c %*% b2[[2]]+ a2
    eta3 = X2.c %*% b2[[3]]+ a2
    
    ME2.class[k, 1] = mean(  - plogis(eta1) * (1-plogis(eta1)) * b2[[1]][k] )
    ME2.class[k, 4] = mean(  (1 - plogis(eta2)) * plogis(eta1) * (1-plogis(eta1)) *  b2[[1]][k] +
                               plogis(eta1) * (- plogis(eta2) * (1-plogis(eta2))) * b2[[2]][k]  ) 
    ME2.class[k, 7] = mean( (1 - plogis(eta3)) * plogis(eta2) * (plogis(eta1) * (1-plogis(eta1))) * b2[[1]][k] +
                              (1 - plogis(eta3)) * plogis(eta1) * (plogis(eta2) * (1-plogis(eta2))) * b2[[2]][k] + 
                              plogis(eta1) * plogis(eta2) * (-plogis(eta3) * (1-plogis(eta3))) * b2[[3]][k] )
    ME2.class[k, 10] = mean( plogis(eta3) * plogis(eta2) * plogis(eta1) * (1-plogis(eta1)) * b2[[1]][k] +
                               plogis(eta3) * plogis(eta1) * plogis(eta2) * (1-plogis(eta2)) * b2[[2]][k] + 
                               plogis(eta1) * plogis(eta2) * plogis(eta3) * (1-plogis(eta3)) * b2[[3]][k] ) 
    
  }
}