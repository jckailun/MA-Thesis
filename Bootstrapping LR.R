##### BS LR test #####
### Cohort 1 ###

g1 = 2; g2 = 3; n1 = nrow(X1); n2 = nrow(X2)
est11.bs = NULL; LR11.bs = NULL; LR01.bs = NULL; est01.bs = NULL; LR01.null.bs = NULL; LR11.null.bs = NULL
set.seed(1000)

for (b in 1:B){
  
  idx.bs = sample(1:n1, n1, replace = TRUE)
  dummy.bs = dummy
  dummy.bs[[1]] = dummy.bs[[1]][idx.bs, ]
  
  # case 1
  theta.01 = rep(0, 49)
  X1.bs = X1.c[idx.bs, ]
  out = brd(theta.01, list(X1.bs), k=16, cohort = 1, dummy.bs)
  est01.bs = rbind(est01.bs, out$zero)
  LR01.bs = c(LR01.bs, ll(out$zero, k=16, list(X1.bs), cohort = 1, dummy.bs))
  
  out.null = brd.lr(0, cohort = 1 , dummy = dummy.bs)
  LR01.null.bs = c(LR01.null.bs, out.null$loglik[length(out.null$loglik)])

  # case 2
  theta.11 = rep(0, 50)
  gower_dist_11 = daisy(X1.bs, metric = "gower", type = list(symm = c(11, 13:16)))
  diss11 = as.matrix(gower_dist_11)
  #c1 = Cluster_Medoids(diss11, clusters = g1, swap_phase = TRUE, verbose = TRUE, fuzzy = TRUE)
  c1 = pam(diss11, k = 2, diss = TRUE)
  
  cl1 = c1$clustering
  
  idx.g1.d1 = NULL; idx.g2.d1 = NULL
  for (i in 1:n1){
    idx.g1.d1 = c(idx.g1.d1, ifelse(cl1[i]==1, i, 0))
    idx.g2.d1 = c(idx.g2.d1, ifelse(cl1[i]==2, i, 0))
  }
  idx.g1.d1 = idx.g1.d1[idx.g1.d1!=0]
  idx.g2.d1 = idx.g2.d1[idx.g2.d1!=0]
  
  g1.d1 = NULL
  for (i in idx.g1.d1){
    g1.d1 = rbind(g1.d1, X1.bs[i,])
  }
  
  g2.d1 = NULL
  for (i in idx.g2.d1){
    g2.d1 = rbind(g2.d1, X1.bs[i,])
  }
  
  G1.bs = list(g1.d1, g2.d1)
  idx.d1 = list(idx.g1.d1, idx.g2.d1)
  out = brd(theta.11, G1.bs, k=16, cohort = 1, dummy.bs)
  est11.bs = rbind(est11.bs, out$zero)
  LR11.bs = c(LR11.bs, ll(out$zero, G1.bs, k=16, cohort = 1, dummy.bs))
  
  out.null = brd.lr(c(0,0), cohort = 1 , dummy = dummy.bs)
  LR11.null.bs = c(LR11.null.bs, out.null$loglik[length(out.null$loglik)])
  
}

LR01.dist = -2*(LR01.null.bs - LR01.bs) - LR01
plot(density(LR01.dist), main = "Bootstrap distribution of LR", sub = 'Case 1, Cohort 1')
LR01.ts = LR01
LR01.pval = mean(LR01.dist >= LR01.ts) ### nvm

LR11.dist = -2*(LR11.null.bs - LR11.bs) - LR11
plot(density(LR11.dist), main = "Bootstrap distribution of LR", sub = 'Case 2, Cohort 1')
LR11.ts = LR11
LR11.pval = mean(LR11.dist >= LR11.ts) ### nvm

cbc.bs.cohort1 = -2*(LR01.bs - LR11.bs) # - ts.cbc.cohort1 ?
plot(density(cbc.bs.cohort1), main = "Bootstrap distribution of LR", sub = 'Case 2, Cohort 1')
ts.cbc.cohort1 = 2*(est11$loglik[length(est11$loglik)] - est01$loglik[length(est01$loglik)])
pval.cbc.cohort1 = mean(cbc.bs.cohort1 >= ts.cbc.cohort1)
### End Cohort 1 ###

### Case 1 ###
cov01.bs = cov(est01.bs)
mu01.bs = apply(est01.bs, 2, mean)
se01.bs = sqrt(diag(cov01.bs))

pval01.bs = NULL
for (v in 1:49){
  ref = est01.bs[, v] - est01$zero[v]
  summand = mean(abs(ref) >= abs(est01$zero[v]))
  pval01.bs = c(pval01.bs, summand)  ###
}

CI01.bs = NULL
for (v in 1:49){
  l = quantile(est01.bs[,v] - est01$zero[v], 0.025)
  u = quantile(est01.bs[,v] - est01$zero[v], 0.975)
  CI01.bs = rbind(CI01.bs, c(est01$zero[v]-u, est01$zero[v]-l))
}

out01.bs = cbind(est01$zero, se01.bs, CI01.bs, pval01.bs)
# write.csv(out01.bs, file='out01bs.csv')
out01.format.bs = NULL
for (i in 1:16){ out01.format.bs = rbind(out01.format.bs, rbind(out01.bs[i,], out01.bs[i+16,], out01.bs[i+32,]))}
write.csv(out01.format.bs, file = 'out01_formatted_bs.csv')
### End Case 1 ###

### Case 2 ###
cov11.bs = cov(est11.bs)
mu11.bs = apply(est11.bs, 2, mean)
se11.bs = sqrt(diag(cov11.bs))
pval11.bs = NULL
for (v in 1:50){
  ref = est11.bs[, v] - est11$zero[v]
  summand = mean(abs(ref) >= abs(est11$zero[v])) ###
  pval11.bs = c(pval11.bs, summand)
}

CI11.bs = NULL
for (v in 1:50){
  l = quantile(est11.bs[,v] - est11$zero[v], 0.025)
  u = quantile(est11.bs[,v] - est11$zero[v], 0.975) 
  CI11.bs = rbind(CI11.bs, c(est11$zero[v]-u, est11$zero[v]-l))
}

out11.bs = cbind(est11$zero, se11.bs, CI11.bs, pval11.bs)
out11.format.bs = NULL
for (i in 1:16){ out11.format.bs = rbind(out11.format.bs, rbind(out11.bs[i,], out11.bs[i+16,], out11.bs[i+32,]))}
write.csv(out11.format.bs, file = 'out11_formatted_bs.csv')
### End Case 2 ###



