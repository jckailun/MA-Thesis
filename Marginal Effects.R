##### Marginal Effects #####
### Cohort 1 ###
b1 = est01$zero
b1 = list(b1[1:16], b1[17:32], b1[33:48]) 

ME1 = matrix(0, nc = 9, nr = 16); J = matrix(0, nr=16, nc=3)
idj = list(1:3, 4:6, 7:9)
idv = list(1:16, 17:32, 33:48)
a = est01$zero[49]
for (j in 1:3){
  for (k in 1:16){
        eta = X1.c %*% b1[[j]] + a
        J[k, j] = mean( plogis(eta) * (1-plogis(eta)) * b1[[j]][k])
  }
}
#rec = NULL
#for (j in 1:3){
#rec = c(rec, sqrt((J[,j]) %*% v1[idv[[j]], idv[[j]]] %*% (J[,j]))) }


### Cohort 2 ###
b2 = est02$zero
b2 = list(b2[1:16], b2[17:32], b2[33:48])
J2 = matrix(0, nr=16, nc=3)
idj = list(1:3, 4:6, 7:9)
idv = list(1:16, 17:32, 33:48)
a2 = est02$zero[49]
for (j in 1:3){
  for (k in 1:16){
    eta = X2.c %*% b2[[j]] + a2
    J2[k, j] = mean( plogis(eta) * (1-plogis(eta)) * b2[[j]][k])
  }
}

write.csv(ME1, 'Marginal Effects of Cohort 1.csv')
write.csv(ME2, 'Marginal Effects of Cohort 2.csv')

### Odds ratios ###
oddsratio1 = matrix(0, nr = 5, nc = 9)
for (j in 1:3){
  for (k in c(11, 13:16)){
    dat.temp = X1.c
    
    dat.temp[, k] = 1
    eta1 = dat.temp %*% b1[[j]]
    eta1.up = dat.temp %*% up01[[j]] 
    eta1.low = dat.temp %*% low01[[j]]
    
    dat.temp[, k] = 0
    eta0 = dat.temp %*% b1[[j]] 
    eta0.up = dat.temp %*% up01[[j]] 
    eta0.low = dat.temp %*% low01[[j]]
    
    oddsratio1[match(k, c(11, 13:16)), idj[[j]][1]] = mean( (plogis(eta1) / (1-plogis(eta1))) / (plogis(eta0) / (1-plogis(eta0))) )
    oddsratio1[match(k, c(11, 13:16)), idj[[j]][2]] = mean( (plogis(eta1.low) / (1-plogis(eta1.low))) / (plogis(eta0.low) / (1-plogis(eta0.low))) )
    oddsratio1[match(k, c(11, 13:16)), idj[[j]][3]] = mean( (plogis(eta1.up) / (1-plogis(eta1.up))) / (plogis(eta0.up) / (1-plogis(eta0.up))) )
    
   }
}


col = c(11, 13:16)
k = permutations(n = 5, r = 2, v = col, repeats.allowed = FALSE)
k = k[c(1:4, 6:8, 11:12, 16), ]
oddsratio.inter1 = matrix(0, nr = nrow(k), nc = 3)

for (j in 1:3){
  for (comb in 1:nrow(k)){
    
    dat.temp = X1.c
    dat.temp[, k[comb,]] = 1
    eta1 = dat.temp %*% b1[[j]]
    dat.temp[, k] = 0
    eta0 = dat.temp %*% b1[[j]] 
    
    oddsratio.inter1[comb, j] = mean( (plogis(eta1) / (1-plogis(eta1))) / (plogis(eta0) / (1-plogis(eta0))) )
  }
}
oddsratio.inter1 = cbind(k, oddsratio.inter1)
nom = c('male', 'nuclear', 'urban', 'black', 'hispanic')

oddsratio2 = matrix(0, nr = 5, nc = 9)
for (j in 1:3){
  for (k in c(11, 13:16)){
    dat.temp = X2.c
    
    dat.temp[, k] = 1
    eta1 = dat.temp %*% b2[[j]]
    eta1.up = dat.temp %*% up02[[j]] 
    eta1.low = dat.temp %*% low02[[j]]
    
    dat.temp[, k] = 0
    eta0 = dat.temp %*% b2[[j]] 
    eta0.up = dat.temp %*% up02[[j]] 
    eta0.low = dat.temp %*% low02[[j]]
    
    oddsratio2[match(k, c(11, 13:16)), idj[[j]][1]] = mean( (plogis(eta1) / (1-plogis(eta1))) / (plogis(eta0) / (1-plogis(eta0))) )
    oddsratio2[match(k, c(11, 13:16)), idj[[j]][2]] = mean( (plogis(eta1.low) / (1-plogis(eta1.low))) / (plogis(eta0.low) / (1-plogis(eta0.low))) )
    oddsratio2[match(k, c(11, 13:16)), idj[[j]][3]] = mean( (plogis(eta1.up) / (1-plogis(eta1.up))) / (plogis(eta0.up) / (1-plogis(eta0.up))) )
    
  }
}



write.csv(oddsratio1, 'Odds Ratio of Cohort 1.csv')
write.csv(oddsratio2, 'Odds Ratio of Cohort 2.csv')

